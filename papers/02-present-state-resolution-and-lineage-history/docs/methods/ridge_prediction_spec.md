# Full ridge coefficient and prediction equations

Let X be molecular coordinates after column-wise training-derived centering/scaling (population standard deviation; replace zero scale by one). Let Y be the training target coordinates, and let Ybar be their training mean. Fit

`B_alpha = solve(X.T @ X + alpha I, X.T @ (Y - Ybar))`.

The actual held-out prediction, **including its intercept**, is

`Yhat(x) = Ybar + x_standardized @ B_alpha`.

For history, transform the first six nonnegative tree coordinates by log1p and leave edit fraction unchanged. Learn all seven column centers/scales from fit-owned cells only. Let X_H = [X | H] be the concatenated standardized design. Fit

`B_alpha,H = solve(X_H.T @ X_H + alpha I, X_H.T @ (Y - Ybar))`,

then predict

`Yhat_H(x,h) = Ybar + [x_standardized,h_standardized] @ B_alpha,H`.

The intercept is unpenalized; all molecular and history coefficients share the selected alpha penalty. Molecular-only and history models select alpha independently using their mean of two inner-embryo NSEs. Original alpha = {0.01,0.1,1,10}; E1 additionally includes {100,1000,10000}. Predictions are row-L2-normalized for NSE/cosine/coordinate MAE. Training scaling and target mean are reapplied unchanged to validation/test cells.

The null implementation uses the algebraically equivalent seven-column Schur complement of these normal equations, with cached molecular Cholesky factorizations. It refits standardized history and repeats every alpha candidate and nested selection under every reassignment; direct concatenated-design solves provide independent numerical checks. It is not a frozen-parameter permutation approximation.
