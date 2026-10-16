# IMM: the mixing step, derived

Notation: `r` motion models; `M_j(k)` means "model `j` is in force during scan `k`". The model sequence is a Markov chain with `pi_ij = P(M_j(k) | M_i(k-1))`. After scan `k-1` each model's filter holds `(x_i, P_i)` and a probability `mu_i = P(M_i(k-1) | Z^{k-1})`, where `Z^{k-1}` is all data so far.

## Why an approximation is needed

The exact Bayes filter for a jump-Markov system conditions on the whole model history. Its posterior is a Gaussian mixture with `r^k` components after `k` scans, which is infeasible. Generalised pseudo-Bayesian filters keep the mixture bounded by merging components. The IMM merges once per model per scan, which leaves `r` filters per scan. GPB1 merges everything into one; GPB2 keeps `r^2`.

## Step 1: predicted model probabilities

By the law of total probability,

```text
c_j = P(M_j(k) | Z^{k-1}) = sum_i pi_ij * mu_i
```

## Step 2: mixing probabilities

Given that model `j` is in force at `k`, which model was in force at `k-1`? Bayes' rule:

```text
mu_{i|j} = P(M_i(k-1) | M_j(k), Z^{k-1}) = pi_ij * mu_i / c_j
```

These weights sum to one over `i`.

## Step 3: mixed initial condition for model j

Conditioned on `M_j(k)`, the state at `k-1` is distributed as the mixture

```text
p(x | M_j(k), Z^{k-1}) = sum_i mu_{i|j} * N(x; x_i, P_i)
```

Model `j`'s filter needs a single Gaussian to start from. Moment matching picks the Gaussian with the same mean and covariance as the mixture.

Mean:

```text
x0_j = sum_i mu_{i|j} * x_i
```

Covariance, by the law of total variance (mean of the covariances plus covariance of the means):

```text
P0_j = sum_i mu_{i|j} * [ P_i + (x_i - x0_j)(x_i - x0_j)^T ]
```

The second term is the spread of the means. Dropping it makes the filter overconfident exactly when the models disagree, which is when a manoeuvre is starting. It is the term that lets the IMM widen its uncertainty at a model switch.

## Step 4: model-matched prediction and update

Run filter `j` (here an EKF) from `(x0_j, P0_j)`: predict with model `j`'s transition `F_j` and process noise `Q_j`, then update with the measurement. Each filter produces an innovation `y_j` and covariance `S_j`.

## Step 5: model probability update

The likelihood of the measurement under model `j` is the Gaussian density of its innovation:

```text
Lambda_j = N(y_j; 0, S_j) = exp(-0.5 * y_j^T S_j^{-1} y_j) / (2*pi*sqrt(det S_j))
mu_j(k)  = Lambda_j * c_j / sum_l Lambda_l * c_l
```

The implementation works in the log domain and subtracts the maximum before exponentiating, so a large innovation cannot underflow every likelihood to zero.

## Step 6: combined output

```text
x = sum_j mu_j * x_j
P = sum_j mu_j * [ P_j + (x_j - x)(x_j - x)^T ]
```

The combined estimate is only for output. Each filter continues from its own `(x_j, P_j)` and goes through Step 3 again next scan.

## Missed detections

If no measurement is assigned on a scan, the update is skipped: each filter keeps its predicted state and `mu_j = c_j`. Model probabilities then drift toward the stationary distribution of `pi`.

## What this implementation assumes

- Models share the state `[x, y, vx, vy]`, so mixing is a direct weighted sum. Models with different state dimensions need a mapping step.
- The coordinated-turn models use fixed turn rates (`+omega` and `-omega`). With a fixed rate the transition is linear. Estimating the turn rate as part of the state would make it nonlinear and is not implemented.
- The Gaussian moment matching and the EKF linearisation are approximations, so the likelihoods are only approximately correct.

## Reference

Blom, H.A.P. and Bar-Shalom, Y., "The interacting multiple model algorithm for systems with Markovian switching coefficients", IEEE Trans. Automatic Control 33(8), 1988.
