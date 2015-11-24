// File:    test/unit_test/TestFPExceptionControl.cpp
// Purpose: test an uniform mechanism to enable/disable Floating Point exceptions

#include <cstdio>   // printf
#include <cstdlib>  // exit()
#include "base/fp_exception_glibc_extension.h"

int main (int argc, char **argv) {
    double s;
    struct sigaction act;

    act.sa_sigaction = (void(*))fhdl;
    sigemptyset (&act.sa_mask);
    act.sa_flags = SA_SIGINFO;

    printf ("Old divByZero exception: 0x%08X\n", feenableexcept (FE_DIVBYZERO));
    printf ("Old invalid exception:   0x%08X\n", feenableexcept (FE_INVALID));
    printf ("New fp exception:        0x%08X\n", fegetexcept ());

    // set handler
    if (sigaction(SIGFPE, &act, (struct sigaction *)0) != 0)
    {
        perror("Yikes");
        exit(-1);
    }

    s = 1.0 / 0.0;  // FE_DIVBYZERO
    s = 0.0 / 0.0;  // FE_INVALID
    return 0;
}
