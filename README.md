# Parametric Curve Parameter Estimation

## Problem Statement

The assignment is to solve a parametric curve represented by the following parametric equations:

x(t) = t·cos(θ) − e^(M|t|)·sin(0.3t)·sin(θ) + X
y(t) = 42 + t·sin(θ) + e^(M|t|)·sin(0.3t)·cos(θ)

for the domain 6 < t < 60. The exact shape and position of this curve can be determined from three unknown parameters:

- **θ** (theta), an angle constrained to 0° < θ < 50°
- **M** is a small exponential growth factor constrained between −0.05 < M < 0.05
- **X** is the horizontal offset that is limited to the range 0 to 100.

A file, called `xy_data.csv`, contained a list of points (x, y) known to be on this curve for the values of t in the given range. Most importantly, the t values corresponding to the points were not provided in the CSV — only the resulting (x, y) values were provided. The problem is to write a program that uses only these points to recover the values of θ, M, and X as closely as possible, not the program itself, but the values of the three parameters (expressed as a LaTeX equation and/or a Desmos link).

## Approach

### Step 1: Recognizing the underlying structure

It was useful to first try to picture what the equations geometrically represent before trying to do any numerical fitting. The system can be written in a matrix form and the structure became much clear and simpler:

$$
\begin{bmatrix} 
x - X \\ 
y - 42 
\end{bmatrix} 
= 
\begin{bmatrix} 
\cos(\theta) & -\sin(\theta) \\ 
\sin(\theta) & \cos(\theta) 
\end{bmatrix} 
\begin{bmatrix} 
t \\ 
e^{M|t|} \cdot \sin(0.3t) 
\end{bmatrix}
$$

This is just a 2-D rotation matrix applied to a more basic curve `(t, e^(M|t|)·sin(0.3t))`, plus a translation of `(X, 42)`. That is, once the rotation and the translation are removed, what remains is a simple exponentially-modulated sine wave, plotted against t on one axis, and t itself on the other. It's crucial to realize this early, since it wasn't just a problem of fitting an arbitrary 2D curve, but more along the lines of finding the rotation angle and offset for the data to re-normalize it back into this known, simple functional form.


### Step 2: First attempt — nearest-neighbor grid matching

The first working approach didn't directly use the structure above, but considered the problem as a generic point-cloud matching problem:

1. For a given candidate set of parameters (θ, M, X), compute a dense set of points along a candidate curve by calculating x(t) and y(t) at a fixed set of points t (e.g., 2000 evenly spaced t's between 6 and 60).
2. At each point in the actual data (in xy_data.csv), determine the point on this candidate curve that is closest to it (in the sense of the sum of the absolute differences in x and y — the L1 distance).
3. Sum all these minimum distances to obtain a total "loss" on this candidate parameter set.
4. Find the values of (θ, M, X) by minimizing this total loss using `scipy.optimize.minimize`, while keeping the values of **θ** between 0 and 50°, and keeping the values of **M** between −0.05 and 0.05, and keeping the values of **X** between 0 and 100.

This was a reasonable method that was easily understood by a lot of students, but there were two practical drawbacks. First, the candidate curve will only be computed at a finite number of points on the grid, and the "closest point" returned will only be as precise as the precision of the grid: the more coarse the grid, the more approximation error will be embedded in all the distance calculations. Second, the `min()` function applied to a discrete array is not smooth, so that small changes in the parameters can lead to big jumps on which point is "closest," causing gradient based optimizers to perform in an unpredictable manner and possibly converging to a "good enough" solution and not the minimum.

### Step 3: Refined approach — analytic residual fitting

The geometric knowledge gained in Step 1 lead to a far more exact approach: Simply apply the inverse rotation and translation to each individual data point, and determine by calculation if it is a solution of the known underlying function.

Mathematically, if (x, y) is any real data point, and (θ, M, X) is any candidate, the data point can be rotated backwards and de-translated.

$$
\begin{aligned}
u &= (x - X) \cdot \cos(\theta) + (y - 42) \cdot \sin(\theta) \\
v &= -(x - X) \cdot \sin(\theta) + (y - 42) \cdot \cos(\theta)
\end{aligned}
$$

If (θ, M, X) are correctly calculated, then the value of`u` should be equal to the point's original (unknown) t-value, and `v` should be equal to `e^(M|u|)·sin(0.3u)`. This provides a “residual” for each data point:

$$
\text{residual} = v - e^{M|u|} \cdot \sin(0.3u)
$$


These residuals were reduced with the use of `scipy.optimize.least_squares`, which uses a trust-region/Levenberg-Marquardt algorithm that is appropriate for smooth nonlinear least-squares. This approach is superior to Step 2 for three reasons: it avoids the use of a grid or discretization, thus avoiding any approximation error due to finite resolution of the grid or finite radius of the ball; the residual function is smooth and differentiable; and the optimizer can be pushed to very tight tolerances, 1e-15 for `xtol`, `ftol`, and `gtol`, and local optima possibilities remain, as the fit was repeated from a range of different starting guesses within the legitimate ranges of the parameters.

### Step 4: Validation

Two independent checks were used to confirm the fitted parameters were correct:

**L1 distance validation.** An initial attempt was to match data points to a curve that was predicted by index (i.e. assume that data row `i` is the i-th value in a uniformly-spaced t-grid). On inspection this gave an unexpectedly large L1 distance (mean ≈ 25), which was found to be due to the x-values in `xy_data.csv` being not monotonically increasing from row to row. The switching approach (using the nearest point on a finely sampled version of the fitted curve instead of assuming any order for the rows) was much more meaningful – with a mean L1 distance of ~0.0041 and a maximum of ~0.0124 at all points, which is truly negligible on the scale of this data (coordinates in the range of ~50-100 units).




**Plotting.**  A fitted curve was plotted with Matplotlib and a scatter plot of the actual data points from `xy_data.csv` was also plotted. The fitted curve line visually passes through every data point, with no visible deviation from the points, further reinforcing the above numerical L1 result.

## Final Results

The fitted parameters converged to:

- **θ ≈ 0.5235983031753396 radians** (≈ 30°)
- **M ≈ 0.029999996873059678**
- **X ≈ 54.99999821279947**

The near-round values (30°, 0.03, 55) resulting from an unconstrained numerical fit also indicate a strong possibility that these are the true, exact parameters with which the original dataset was produced; the slight deviation from perfectly round numbers is due to rounding and/or the tolerance of the optimizer, and not to any real error in the numerical fit.

## Final Equation

\left(t*\cos(0.5235983031753396)-e^{0.029999996873059678\left|t\right|}\cdot\sin(0.3t)\sin(0.5235983031753396)+54.99999821279947,42+t*\sin(0.5235983031753396)+e^{0.029999996873059678\left|t\right|}\cdot\sin(0.3t)\cos(0.5235983031753396)\right)

Domain: 6 ≤ t ≤ 60

**Desmos link:** https://www.desmos.com/calculator/me0yi3twh7

## Tools & References
- [SciPy `least_squares` documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html) — used for parameter optimization (Levenberg-Marquardt / trust-region reflective algorithm)
- [NumPy](https://numpy.org/) — array computation and vectorized math
- [pandas](https://pandas.pydata.org/) — CSV loading and data handling
- [Matplotlib](https://matplotlib.org/) — visualization of fitted curve vs. actual data
- [Desmos](https://www.desmos.com/) — interactive curve graphing and equation verification