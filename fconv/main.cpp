#include "fkrelu.h"
#include "octahedron.h"
#include "pdd.h"
#include "split_in_quadrants.h"
#include "utils.h"
#include <execinfo.h>
#include <iostream>
#include <map>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

using namespace std;

using namespace Eigen;

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

void run_fkrelu(int k, int i) {
  string path =
      "octahedron_hrep/k" + to_string(k) + "/" + to_string(i) + ".txt";
  cout << "working with path: " << path << endl;
  MatrixXd hrep = read_matrix(path);
  Timer t;
  MatrixXd res = fkrelu(hrep);
  int micros = t.micros();
  cout << "result contains " << res.rows() << " rows and took " << micros / 1000
       << " ms" << endl;
}

// Looks like there are a bunch of segfault errors.
int main() {
  signal(SIGSEGV, handler);
  dd_set_global_constants();

  cout << "before cycle" << endl;

  run_fkrelu(2, 1);

  for (int k = 2; k <= 4; k++) {
    for (int i = 1; i <= 4; i++) {
      run_fkrelu(k, i);
    }
  }

  dd_free_global_constants();
  return 0;
}
