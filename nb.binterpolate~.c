/*
    nb.binterpolate~:
    Signal processing external written for Max version 7.1.0
 
    nb.binterpolate~ takes an FFT signal as input and interpolates between the values of each bin measured at different points in time.
    
    Naithan Bosse - 2017
    https://naithan.com
*/
 
#include "ext.h"
#include "ext_obex.h"
#include "z_dsp.h"
#include "r_pfft.h"

#define DEFAULT_FFT_SIZE 4096
#define DEFAULT_LENGTH 10
#define MAX_LENGTH 30
#define MIN_LENGTH 0
#define DEFAULT_VARIANCE 2
#define MAX_VARIANCE 15
#define MIN_VARIANCE 0
#define MIN_INTERP_FRAMES 1

typedef struct _interp {
	t_pxobject	ob;             // The object "base class"
	double		fftSize;
    int         sampleRate;

    t_atom*     currMag;        // Contains the current magnitude/real values for each FFT bin (used while interpolating to the target magnitudes)
    t_atom*     currPhase;      // Contains the current phase/imaginary values for each FFT bin (used interpolating to the target phases)
    t_atom*     targetMag;      // Target list of magnitudes for the interpolation
    t_atom*     targetPhase;    // Target list of phases for the interpolation
    t_atom*     incMag;         // Amount to increment each magnitude per frame
    t_atom*     incPhase;       // Amount to increment each phase per frame
    t_atom*     totalFrames;    // Total number of frames used for the interpolation
    t_atom*     frameCount;     // The current frame of the interpolation (frameCount/totalFrames * 100 = interpolation %)
    t_atom*     updateTarget;   // For each FFT bin, set to updateTarget to 1 if interpolation is complete and we need a new target, 0 otherwise
    
    float       interpLengthSecs;       // Base number of seconds to spend interpolating between the start and goal FFT snapshots
    float       interpVarianceSecs;     // Maximum allowable random variance to add or subtract to the base interpolation length (in seconds)
    int         interpLengthFrames;     // (interpLengthSecs * sampleRate) / fftSize
    int         interpVarianceFrames;   // (interpVarianceSecs * sampleRate) / fftSize
    int         interpMin;              // Max((interpLengthFrames - interpVarianceFrames), 1)
    int         interpMax;              // interpLengthFrames + interpVarianceFrames
} t_interp;


// Method prototypes
void *interp_new(t_symbol *s, long argc, t_atom *argv);
void interp_free(t_interp *x);
void interp_assist(t_interp *x, void *b, long m, long a, char *s);
void interp_bang(t_interp *x);
void interp_float(t_interp *x, double f);
void interp_int(t_interp *x, long n);
void interp_dsp64(t_interp *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);
void interp_perform64(t_interp *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);

// Helper functions
void updateTarget(t_interp *x, t_double mag, t_double phase, t_double fftBin);
double getFFTSize(t_interp *x);
void setInterpolationTime(t_interp *x, float interpLengthSecs, float interpVarianceSecs);
int secondsToFrames(float seconds, int sampleRate, int fftSize);
float frand(float min, float max);
int irand(int min, int max);

// Global class pointer variable
static t_class *interp_class = NULL;

//***********************************************************************************************
// Max class methods
//***********************************************************************************************

/**
 * Main function for the external
 * Registers the class and adds the appropriate class methods
 */
void ext_main(void *r) {
	t_class *c = class_new("nb.binterpolate~", (method)interp_new, (method)interp_free, (long)sizeof(t_interp), 0L, A_GIMME, 0);

    class_addmethod(c, (method)interp_bang,     "bang",                 0);
    class_addmethod(c, (method)interp_int,      "int",      A_LONG,     0);
    class_addmethod(c, (method)interp_float,    "float",    A_FLOAT,    0);
	class_addmethod(c, (method)interp_dsp64,	"dsp64",	A_CANT,     0);
	class_addmethod(c, (method)interp_assist,	"assist",	A_CANT,     0);

	class_dspinit(c);
	class_register(CLASS_BOX, c);
	interp_class = c;
}


/**
 * Constructor
 */
void *interp_new(t_symbol *s, long argc, t_atom *argv) {
	t_interp *x = (t_interp *)object_alloc(interp_class);
	if (x) {
		dsp_setup((t_pxobject *)x, 3);	// MSP inlets: argument 2 is the # of inlets
        x->ob.z_misc = Z_NO_INPLACE;
        x->sampleRate = sys_getsr();
        x->fftSize = getFFTSize(x);
        srandom(time(NULL)); // Seed random numbers with the time the object is created

        // Create outlets
        outlet_new(x, "signal");
        outlet_new(x, "signal");
        
        // Set interpolation length and variance using arguments if available
        float interpLength = (argc > 0) ? atom_getfloat(argv) : DEFAULT_LENGTH;
        float interpVariance = (argc > 1) ? atom_getfloat(argv+1) : DEFAULT_VARIANCE;
        setInterpolationTime(x, interpLength, interpVariance);
        
		// Allocate memory
        x->currMag      = (t_atom*)sysmem_newptrclear(sizeof(t_atom) * x->fftSize);
        x->currPhase    = (t_atom*)sysmem_newptrclear(sizeof(t_atom) * x->fftSize);
        x->targetMag    = (t_atom*)sysmem_newptrclear(sizeof(t_atom) * x->fftSize);
        x->targetPhase  = (t_atom*)sysmem_newptrclear(sizeof(t_atom) * x->fftSize);
        x->incMag       = (t_atom*)sysmem_newptrclear(sizeof(t_atom) * x->fftSize);
        x->incPhase     = (t_atom*)sysmem_newptrclear(sizeof(t_atom) * x->fftSize);
        x->totalFrames  = (t_atom*)sysmem_newptrclear(sizeof(t_atom) * x->fftSize);
        x->frameCount   = (t_atom*)sysmem_newptrclear(sizeof(t_atom) * x->fftSize);
        x->updateTarget = (t_atom*)sysmem_newptrclear(sizeof(t_atom) * x->fftSize);
	}
	return (x);
}

/**
 * Destructor
 */
void interp_free(t_interp *x) {
    dsp_free((t_pxobject*)x);
    sysmem_freeptr(x->currMag);
    sysmem_freeptr(x->currPhase);
    sysmem_freeptr(x->targetMag);
    sysmem_freeptr(x->targetPhase);
    sysmem_freeptr(x->incMag);
    sysmem_freeptr(x->incPhase);
}

/**
 * Show assist messages when mouse hovers over each inlet/outlet.
 */
void interp_assist(t_interp *x, void *b, long m, long a, char *s) {
	if (m == ASSIST_INLET) {
        if (a == 0) {
            sprintf(s, "(Signal) FFT magnitude or real component\n(Float) Interpolation length in seconds");
        } else if(a == 1) {
            sprintf(s, "(Signal) FFT phase or imaginary component\n(Float) Interpolation variance in seconds");
        } else if(a == 2) {
            sprintf(s, "(Signal) FFT index");
        }
	} else if (m == ASSIST_OUTLET) {
        if (a == 0) {
            sprintf(s, "(Signal) FFT magnitude or real component");
        } else if (a == 1) {
            sprintf(s, "(Signal) FFT phase or imaginary component");
        }
	}
}

/**
 * Handle bang (just prints to the Max console at the moment)
 * @param x pointer to the object struct
 */
void interp_bang(t_interp *x) {
    post("nb.binterpolate~ was written by Naithan Bosse in 2017");
}

/**
 * Handle float input
 * @param x pointer to the object struct
 * @param f the interpolation length or variance
 * If the float is sent in input 0, set the interpolation length. If the float is sent in input 1, set the interpolation variance.
 */
void interp_float(t_interp *x, double f) {
    long inlet = proxy_getinlet((t_object *)x);
    if (inlet == 0) {
        // Inlet 0: Use input to set interpolation length
        setInterpolationTime(x, f, x->interpVarianceSecs);
    } else if (inlet == 1) {
        // Inlet 1: Use input to set the random variance added to the interpolation length
        setInterpolationTime(x, x->interpLengthSecs, f);
    }
}

/**
 * Handle integer input
 * @param x pointer to the object struct
 * @param f the interpolation length or variance
 */
void interp_int(t_interp *x, long n) {
    interp_float(x,(double)n);
}

//***********************************************************************************************
// Helper functions
//***********************************************************************************************
/**
 * Random number helper function
 * @return a random float between min and max
 */
float frand(float min, float max) {
    float scale = (float)random()/(float)RAND_MAX;
    return min + (scale * (max-min));
}

/**
 * Random number helper function
 * @return a random int between min and max
 */
int irand(int min, int max) {
    float scale = (float)random()/(float)RAND_MAX;
    return min + ((float)scale * (float)(max-min));
}

/**
 * Get the fft size from a pfft~ object containing nb.binterpolate~
 * Use default FFT size if nb.binterpolate~ is not inside a pfft~ object
 */
double getFFTSize(t_interp *x) {
    t_pfftpub *pfft = (t_pfftpub*)gensym("__pfft~__")->s_thing;
    if (pfft)
        return pfft->x_fftsize;
    else
        return DEFAULT_FFT_SIZE;
}

/**
 * @return number of frames
 */
int secondsToFrames(float seconds, int sampleRate, int fftSize) {
    return (int)((seconds * sampleRate) / fftSize);
}

/**
 * Set the minimum and maximum interpolation times that will be randomly chosen in the perform method.
 */
void setInterpolationTime(t_interp *x, float interpLengthSecs, float interpVarianceSecs) {
    // Set the base interpolation length in both seconds and frames
    x->interpLengthSecs = CLAMP(interpLengthSecs, MIN_LENGTH, MAX_LENGTH);
    x->interpLengthFrames = secondsToFrames(x->interpLengthSecs, x->sampleRate, x->fftSize);
    x->interpLengthFrames = (x->interpLengthFrames < MIN_INTERP_FRAMES) ? MIN_INTERP_FRAMES : x->interpLengthFrames;
    
    // Set the amount of random variance in both seconds and frames
    x->interpVarianceSecs = CLAMP(interpVarianceSecs, MIN_VARIANCE, MAX_VARIANCE);
    x->interpVarianceFrames = secondsToFrames(x->interpVarianceSecs, x->sampleRate, x->fftSize);
    
    // Set the min/max frame values
    double minVar = x->interpLengthFrames-x->interpVarianceFrames;
    x->interpMin = (minVar <= 0) ? 1 : minVar;
    x->interpMax = x->interpLengthFrames+x->interpVarianceFrames;
}

/**
 * Update the target value for the interpolation using the given magnitude and phase values
 * @param x pointer to the object struct
 * @param mag the new magnitude value
 * @param phase the new phase value
 * @param fftBin the fftBin index to update
 */
void updateTarget(t_interp *x, t_double mag, t_double phase, t_double fftBin) {
    // Set interpolation target to current signal value for the current fft bin
    int bin = CLAMP(fftBin, 0, x->fftSize-1);
    atom_setfloat(x->targetMag+bin, mag);
    atom_setfloat(x->targetPhase+bin, phase);
    
    // Calculate how much to increment the current bin each frame
    atom_setlong(x->totalFrames+bin, irand(x->interpMin, x->interpMax));
    double magDelta = (mag-atom_getfloat(x->currMag+bin));
    double phaseDelta = (phase-atom_getfloat(x->currPhase+bin));
    magDelta = magDelta/atom_getlong(x->totalFrames+bin);
    phaseDelta = phaseDelta/atom_getlong(x->totalFrames+bin);
    atom_setfloat(x->incMag+bin, magDelta);
    atom_setfloat(x->incPhase+bin, phaseDelta);
    
    // Reset updateTarget flag and counter
    atom_setlong(x->updateTarget+bin, 0);
    atom_setlong(x->frameCount+bin, 0);
}

//***********************************************************************************************
// DSP
//***********************************************************************************************
/**
 * Setup DSP, called when audio is activated
 * Registers the 64-bit perform method in the signal chain in MSP
 */
void interp_dsp64(t_interp *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags) {
    x->sampleRate = samplerate; // Update the sample rate in case it has changed since the object was created
    
    // Set updateTarget to true when audio is started so that we get a new interpolation target.
    for (int i = 0; i<x->fftSize; i++) {
        atom_setlong(x->updateTarget+i, 1);
    }

    object_method(dsp64, gensym("dsp_add64"), x, interp_perform64, 0, NULL);
}

/**
 * 64-bit audio perform method
 */
void interp_perform64(t_interp *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam) {
    // Input signal vectors
	t_double *in_mag = ins[0];      // Either the real component of the FFT or the magnitude depending on if the input is in cartesian or polar coordinates
    t_double *in_phase = ins[1];    // Either the imaginary component of the FFT or the phase depending on the coordinate system
    t_double *in_index = ins[2];    // The FFT bin index
    
    // Output signal vectors
    t_double *out_mag = outs[0];	// Left outlet - magnitude/real
    t_double *out_phase = outs[1];  // Right outlet - phase/imaginary
    
    long n = sampleframes;          // Signal vector size

    // Update target if interpolation is complete.
    for (int i=0; i<sampleframes; i++) {
        if (x->updateTarget+i != NULL && atom_getlong(x->updateTarget+i)) {
            // Target reached - Set the old target value as the new starting point for interpolation
            atom_setfloat(x->currMag+i, atom_getfloat(x->targetMag+i));
            atom_setfloat(x->currPhase+i, atom_getfloat(x->targetPhase+i));
            
            // Update target with new inputs
            updateTarget(x, in_mag[i], in_phase[i], in_index[i]);
        }
    }
    
    // Increment each bin by the appropriate amount then write its value to the corresponding output channel
    int k = 0;
    while (n--) {
        // Get the FFT bin index and CLAMP it between 0 and x->fftSize to avoid a segfault if x->fftSize doesn't match the outer fft size.
        // (Note: This won't occur if the object is used inside a pfft~ object as intended.)
        int bin = CLAMP((int)(in_index[k]), 0, x->fftSize-1);
        k++;
        
        // Get the increment amount for the current bin and add this amount to the current value for the bin
        double incM = atom_getfloat(x->incMag+bin);
        double incP = atom_getfloat(x->incPhase+bin);
        atom_setfloat(x->currMag+bin, atom_getfloat(x->currMag+bin)+incM);
        atom_setfloat(x->currPhase+bin, atom_getfloat(x->currPhase+bin)+incP);
        
        // Write the updated value to the output channels
        *out_mag++ = atom_getfloat(x->currMag+bin);
        *out_phase++ = atom_getfloat(x->currPhase+bin);
        
        // Increment frameCount and set the updateTarget flag to true if the current bin has reached its target
        int framePos = atom_getlong(x->frameCount+bin)+1;
        atom_setlong(x->frameCount+bin, framePos);
        if (framePos >= atom_getlong(x->totalFrames+bin)) {
            atom_setlong(x->updateTarget+bin, 1);
            atom_setlong(x->frameCount+bin, 0);
        }
    }
}
