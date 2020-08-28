#define _USE_MATH_DEFINES

#include "dsmc.h"
#include "log.h"
#include "Timer.h"
#include "argHandler.h"

#include "potential.h"
#include "potentials/AxiconBox.h"
#include "potentials/Gravity.h"
#include "potentials/PowerLaw.h"
#include "potentials/Gaussian.h"
#include "potentials/Harmonic.h"
#include "potentials/CrossedBeam.h"

int main(int argc, char const *argv[])
{
    Potentials::Register("axicon_box", &AxiconBox::Create);
    Potentials::Register("gravity", &Gravity::Create);
    Potentials::Register("power_law", &PowerLaw::Create);
    Potentials::Register("crossed_beam", &CrossedBeam::Create);
    
    args::handle(argc, argv, {
        {"-ow", DSMC::overwrite},
        {"-animate", DSMC::animate},
        {"-seed", DSMC::randomize},
        {"-force", DSMC::continueForHighPSD},
        {"-tas", DSMC::use_subCells}
    });
    
    Timer t = Timer("The whole script (initialisation + simulation)");
    
    DSMC sim;
    sim.run();
}
