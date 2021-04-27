#ifndef BEESWARM_BEESWARM_H
#define BEESWARM_BEESWARM_H

void calculateSwarm(double *x, int *n, int *priority,
    int *order, int *side, int *placed, double *workspace, double *y);

void calculateCompactSwarm(double *x, int *n, int *side,
    int *placed, double *workspace, double *y);

#endif /* BEESWARM_BEESWARM_H */
