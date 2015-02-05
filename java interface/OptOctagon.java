package apron;

/**
 * Manager factory for the OptOctagon abstract domain.
 */
public class OptOctagon 
    extends Manager
{

    // Internals
    ////////////

    private native void init();

    private static native void class_init();

    static { System.loadLibrary("japron"); class_init(); }


    // Constructors
    ///////////////
    
    /**
     * Creates a new manager to create and manipulate OptOctagons.
     */
    public OptOctagon()
    {
        init(); 
    }

}
