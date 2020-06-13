#
#
#  This source file is part of ELINA (ETH LIbrary for Numerical Analysis).
#  ELINA is Copyright © 2019 Department of Computer Science, ETH Zurich
#  This software is distributed under GNU Lesser General Public License Version 3.0.
#  For more information, see the ELINA project website at:
#  http://elina.ethz.ch
#
#  THE SOFTWARE IS PROVIDED "AS-IS" WITHOUT ANY WARRANTY OF ANY KIND, EITHER
#  EXPRESS, IMPLIED OR STATUTORY, INCLUDING BUT NOT LIMITED TO ANY WARRANTY
#  THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS OR BE ERROR-FREE AND ANY
#  IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
#  TITLE, OR NON-INFRINGEMENT.  IN NO EVENT SHALL ETH ZURICH BE LIABLE FOR ANY     
#  DAMAGES, INCLUDING BUT NOT LIMITED TO DIRECT, INDIRECT,
#  SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN
#  ANY WAY CONNECTED WITH THIS SOFTWARE (WHETHER OR NOT BASED UPON WARRANTY,
#  CONTRACT, TORT OR OTHERWISE).
#
#

import sys
import math
import copy
import torch
import numpy
import random
import torch.nn as nn
import torch.nn.functional as F
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score
from datetime import datetime


use_cuda = 'cuda' if torch.cuda.is_available() else 'cpu'
device = torch.device(use_cuda)


class GraphConvolution(nn.Module):

    def __init__(self, in_features, out_features, bias=True):
        super(GraphConvolution, self).__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.weight = nn.Parameter(torch.FloatTensor(in_features, out_features))
        if bias:
            self.bias = nn.Parameter(torch.FloatTensor(out_features))
        else:
            self.register_parameter('bias', None)
        self.reset_parameters()

    def reset_parameters(self):
        stdv = 1. / math.sqrt(self.weight.size(1))
        self.weight.data.uniform_(-stdv, stdv)
        if self.bias is not None:
            self.bias.data.uniform_(-stdv, stdv)

    def forward(self, input, adj):
        support = torch.mm(input, self.weight)
        output = torch.spmm(adj, support)
        if self.bias is not None:
            return output + self.bias
        else:
            return output

    def __repr__(self):
        return self.__class__.__name__ + ' (' \
               + str(self.in_features) + ' -> ' \
               + str(self.out_features) + ')'


def get_big_graph(x, y, res, edges):
    big_x, big_y, big_res, big_edges = [], [], [], []

    off = 0
    for x1, y1, res1, edges1 in zip(x, y, res, edges):
        big_x += x1
        big_y += y1
        for i, res2 in enumerate(res1):
            if res2 == 1:
                big_res.append(i + off)
        for v1, v2, w in edges1:
            big_edges.append([v1 + off, v2 + off, w])

        off += len(x1)

    return big_x, big_y, big_res, big_edges


def balance_ratio(Y, Res):
    pos_num, neg_num = 0, 0
    for y, res in zip(Y, Res):
        for b, c in zip(y, res):
            if c != 1:
                continue

            if b:
                pos_num += 1
            else:
                neg_num += 1

    total_num = pos_num + neg_num
    print('total_num:', total_num)
    print('pos_num:', pos_num)
    print('neg_num:', neg_num)

    return pos_num, neg_num


class PolicyGCN(nn.Module):

    def __init__(self):
        super().__init__()

        self.input_dim = 12
        self.hidden_dim = 128
        self.output_dim = 2

        self.scaler = StandardScaler()

        self.action_net = nn.Sequential(
            nn.Linear(self.hidden_dim, self.hidden_dim),
            nn.ReLU(),
            nn.Linear(self.hidden_dim, self.hidden_dim),
            nn.ReLU(),
            nn.Linear(self.hidden_dim, self.hidden_dim),
            nn.ReLU(),
            nn.Linear(self.hidden_dim, self.output_dim),
        ).to(device)

        self.gc1 = GraphConvolution(self.input_dim, self.hidden_dim).to(device)
        self.gc2 = GraphConvolution(self.hidden_dim, self.hidden_dim).to(device)
        self.gc3 = GraphConvolution(self.hidden_dim, self.hidden_dim).to(device)


    def forward(self, x, edges):
        len_x = len(x)
        x = self.scaler.transform(x)
        x = torch.FloatTensor(x).to(device)
        vertices = []
        weights = []
        for edge in edges:
            vertices.append([edge[0], edge[1]])
            weights.append(edge[2])
        vertices = torch.LongTensor(vertices).to(device)
        weights = torch.FloatTensor(weights).to(device)
        if use_cuda == 'cuda':
            adj = torch.cuda.sparse.FloatTensor(
                vertices.t(),
                weights,
                torch.Size([len_x, len_x])
            )
        else:
            adj = torch.sparse.FloatTensor(
                vertices.t(),
                weights,
                torch.Size([len_x, len_x])
            )

        x = F.relu(self.gc1(x, adj))
        x = F.relu(self.gc2(x, adj))
        x = F.relu(self.gc3(x, adj))
        scores = self.action_net(x)
        return scores


    def do_train(self, x, y, res, edges):
        pos_num, neg_num =  balance_ratio(y, res)
        total_num = pos_num + neg_num

        x_tmp = []
        for a in x:
            x_tmp += a
        random.shuffle(x_tmp)
        self.scaler = self.scaler.fit(x_tmp)

        t = list(zip(x, y, res, edges))
        random.shuffle(t)
        x[:], y[:], res[:], edges[:] = zip(*t)

        epochs = 0
        BATCH_SIZE = 128
        NUM_EPOCHS = 100
        class_weight = None if neg_num == 0 else torch.FloatTensor([neg_num/total_num, pos_num/total_num]).to(device)
        # class_weight = None
        loss_func = nn.CrossEntropyLoss(weight=class_weight)
        optimizer = torch.optim.Adam(self.parameters(), lr=0.0001, weight_decay=1e-5)
        # optimizer = torch.optim.Adam(self.parameters(), lr=0.001, weight_decay=1e-5)
        # scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, lambda epoch: 1 / (epochs // 50 + 1))
        num_batch = (len(x) + BATCH_SIZE - 1) // BATCH_SIZE
        epoch_loss = 0

        best_accuracy = -1
        best_loss = -1
        best_state_dict = None
        for _ in range(NUM_EPOCHS):
            self.train()
            # scheduler.step()

            epoch_loss = 0
            for j in range(num_batch):
                x_batch = x[j*BATCH_SIZE:(j+1)*BATCH_SIZE]
                y_batch = y[j*BATCH_SIZE:(j+1)*BATCH_SIZE]
                res_batch = res[j*BATCH_SIZE:(j+1)*BATCH_SIZE]
                edges_batch = edges[j*BATCH_SIZE:(j+1)*BATCH_SIZE]
                big_x, big_y, big_res, big_edges = get_big_graph(x_batch, y_batch, res_batch, edges_batch)

                scores = self(big_x, big_edges)

                scores = scores[big_res,]
                big_y = torch.LongTensor(big_y).to(device)
                big_res = torch.LongTensor(big_res).to(device)
                big_y = big_y[big_res]
                loss = loss_func(scores, big_y)
                epoch_loss += loss.item()
                optimizer.zero_grad()
                loss.backward()
                optimizer.step()

            score = 0
            for j in range(num_batch):
                x_batch = x[j*BATCH_SIZE:(j+1)*BATCH_SIZE]
                y_batch = y[j*BATCH_SIZE:(j+1)*BATCH_SIZE]
                res_batch = res[j*BATCH_SIZE:(j+1)*BATCH_SIZE]
                edges_batch = edges[j*BATCH_SIZE:(j+1)*BATCH_SIZE]
                big_x, big_y, big_res, big_edges = get_big_graph(x_batch, y_batch, res_batch, edges_batch)

                big_y_pred = self.predict(big_x, big_edges)
                big_y = numpy.array(big_y)
                big_y_pred = numpy.array(big_y_pred)
                big_res = numpy.array(big_res)
                big_y = big_y[big_res]
                big_y_pred = big_y_pred[big_res]
                score += accuracy_score(big_y, big_y_pred)

            epoch_loss = epoch_loss / num_batch
            accuracy = score / num_batch
            if accuracy - best_accuracy > -1e-6:
                if accuracy - best_accuracy > 1e-6:
                    best_accuracy = accuracy
                    best_loss = loss
                    best_state_dict = copy.deepcopy(self.state_dict())
                elif best_loss - loss > 1e-5:
                    best_accuracy = accuracy
                    best_loss = loss
                    best_state_dict = copy.deepcopy(self.state_dict())

            print('[{}] epoch {}, loss {}, accuracy {}'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), epochs, epoch_loss, accuracy))
            sys.stdout.flush()
            epochs += 1

        self.load_state_dict(best_state_dict)


    def predict(self, x, edges):
        self.eval()
        with torch.no_grad():
            scores = self(x, edges)
            _, y = torch.max(scores, 1)
            return y.tolist()


    @staticmethod
    def load(load_path):
        if use_cuda == 'cuda':
            model_file = torch.load(load_path)
        else:
            model_file = torch.load(load_path, map_location='cpu')
        model = PolicyGCN()
        model.load_state_dict(model_file['state_dict'])
        model.scaler = model_file['scaler']

        return model


    def save(self, save_path):
        torch.save({
            'state_dict': self.state_dict(),
            'scaler': self.scaler,
        }, save_path)