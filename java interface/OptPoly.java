/*
	Copyright 2016 Software Reliability Lab, ETH Zurich

	Licensed under the Apache License, Version 2.0 (the "License");
	you may not use this file except in compliance with the License.
	You may obtain a copy of the License at

		http://www.apache.org/licenses/LICENSE-2.0

	Unless required by applicable law or agreed to in writing, software
	distributed under the License is distributed on an "AS IS" BASIS,
	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
	See the License for the specific language governing permissions and
	limitations under the License.
*/

package apron;

/**
 * Manager factory for the convex polyhedron abstract domain.
 */
public class OptPoly 
    extends Manager
{

    // Internals
    ////////////

    private native void init(boolean strict);

    private static native void class_init();

    static { System.loadLibrary("japron"); class_init(); }


    // Constructors
    ///////////////

    /**
     * Creates a new manager to create and manipulate convex polyhedra.
     *
     * <p> If strict is true, then the domain can express strict
     * inequalities, i.e., non-closed convex polyhedra.
     * Otherwise, the domain expresses closed convex polyhedra.
     */
    public OptPoly(boolean strict)
    {
        init(strict); 
    }

}
