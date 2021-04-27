#include "beeswarm.h"

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <R_ext/Utils.h>
#include <R_ext/Visibility.h>

static bool is_y_feasible(double y, double *pre_dxsq_y,
    int start_idx, int end_idx)
{
  for (int i=start_idx; i<end_idx; i++) {
    double pre_dx_sq = pre_dxsq_y[i * 2];
    double pre_y = pre_dxsq_y[i * 2 + 1];
    double y_diff = y - pre_y;
    if (pre_dx_sq + y_diff * y_diff < 0.999)
      return false;
  }
  return true;
}

int compare_y(const void *a, const void *b)
{
  double y0 = ((double *)a)[1];
  double y1 = ((double *)b)[1];
  return (y0 > y1) - (y0 < y1);
}

bool can_use_poty(int j, double poty, double *nearby_dxsq_y, int nearby_xy_len)
{
  int start = j;
  while (start > 0 && nearby_dxsq_y[2*(start-1)+1] > poty - 1) --start;
  int end = j;
  while (end < nearby_xy_len && nearby_dxsq_y[2*end+1] < poty + 1) ++end;
  return is_y_feasible(poty, nearby_dxsq_y, start, end);
}

static void swarm(double *x, int n, int *priority, int *order, int side,
    int *placed, double *workspace, double *y)
{
  for (int iter=0; iter<n; iter++) {
    // place a point

    int i = order[iter];

    R_CheckUserInterrupt();

    // nearby_dxsq_y is a collection of pairs, stored as a flat array.
    // Each pair corresponds to a previously-placed point that is within a
    // distance of 1 horizontally from the point that will be placed.  The
    // pair contains:
    //    y value
    //    square of horizontal distance between this point and the point to be placed.
    // nearby_xy_len is the number of pairs in nearby_dxsq_y
    double *nearby_dxsq_y = workspace;
    int nearby_xy_len = 0;

    int start = i;
    while (start > 0 && x[start-1] > x[i] - 1) --start;
    for (int j=start; j<n; j++) {
      if (priority[j] >= iter) continue;   // skip unplaced points
      if (x[j] >= x[i] + 1) break;
      double x_diff = x[i] - x[j];
      nearby_dxsq_y[nearby_xy_len * 2] = x_diff * x_diff;
      nearby_dxsq_y[nearby_xy_len * 2 + 1] = y[j];
      nearby_xy_len++;
    }

    qsort(nearby_dxsq_y, nearby_xy_len, 2 * sizeof(double), compare_y);

    if (is_y_feasible(0, nearby_dxsq_y, 0, nearby_xy_len)) {
      y[i] = 0;
    } else {
      y[i] = DBL_MAX;
      for (int j=0; j<nearby_xy_len; j++) {
        double dx_sq = nearby_dxsq_y[j * 2];
        double y_ = nearby_dxsq_y[j * 2 + 1];
        double poty_off = sqrt(1 - dx_sq);
        if (side != -1) {
          double poty = y_ + poty_off;
          if (fabs(poty) <= fabs(y[i]) && can_use_poty(j, poty, nearby_dxsq_y, nearby_xy_len))
            y[i] = poty;
        }
        if (side != 1) {
          double poty = y_ - poty_off;
          if (fabs(poty) < fabs(y[i]) && can_use_poty(j, poty, nearby_dxsq_y, nearby_xy_len))
            y[i] = poty;
        }
      }
    }
  }
}

/* Return the index of a unplaced point whose best permitted y value is as
 * close as possible to zero
 */
static int which_min_abs(double *y_best, int *placed, int n)
{
  int i = 0;
  while (placed[i]) ++i;
  double best_val = fabs(y_best[i]);
  int result = i;
  for (++i; i<n; ++i) {
    if (placed[i]) continue;
    double a = fabs(y_best[i]);
    if (a < best_val) {
      best_val = a;
      result = i;
    }
  }
  return result;
}

static void compactSwarm(double *x, int n, int side,
    int *placed, double *workspace, double *y)
{
  double *y_low = workspace;       // largest permitted negative y value
  double *y_high = workspace + n;  // smallest permitted positive y value
  double *y_best = workspace + 2 * n; //  current best permitted y value

  for (int iter=0; iter<n; iter++) {
    R_CheckUserInterrupt();
    int i = which_min_abs(y_best, placed, n);
    double xi = x[i];
    double yi = y_best[i];
    y[i] = yi;
    placed[i] = 1;
    for (int j=0; j<n; j++) {
      if (placed[j]) continue;
      double xdiff = fabs(xi - x[j]);
      if (xdiff >= 1) continue;
      double y_offset = sqrt(1 - xdiff * xdiff);
      double y_hi = fmax(y_high[j], yi + y_offset);
      y_high[j] = y_hi;
      if (side == 0) {
        double y_lo = fmin(y_low[j], yi - y_offset);
        y_low[j] = y_lo;
        y_best[j] = -y_lo < y_hi ? y_lo : y_hi;
      } else {
        y_best[j] = y_hi;
      }
    }
  }

  if (side == -1)
    for (int i=0; i<n; i++)
      y[i] = -y[i];
}

/* Compute a beeswarm layout for the array x.
 *
 * Parameters:
 * x         circle positions on data axis
 * n         length of x
 * priority  1 for first point to be placed etc.
 * order     inverse permutation of `priority`
 * side      -1, 0, or 1
 * placed    which circles have been placed (logical type)
 * workspace an array of doubles for internal use
 * y         (output) circle positions on non-data axis
 */
void attribute_hidden calculateSwarm(double *x, int *n, int *priority,
    int *order, int *side, int *placed, double *workspace, double *y)
{
  swarm(x, *n, priority, order, *side, placed, workspace, y);
}

/* Compute a compact beeswarm layout for the array x.
 *
 * Parameters:
 * x         circle positions on data axis
 * n         length of x
 * side      -1, 0, or 1
 * placed    which circles have been placed (logical type)
 * workspace an array of doubles for internal use
 * y         (output) circle positions on non-data axis
 */
void attribute_hidden calculateCompactSwarm(double *x, int *n, int *side,
    int *placed, double *workspace, double *y)
{
  compactSwarm(x, *n, *side, placed, workspace, y);
}
