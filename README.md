In-progress manuscript, code, and data for the paper tentatively titled "Dealing with Separation in Logistic Regression Models."

Here's the main point: 

When dealing with separation under logistic regression, you had better choose your prior distribution carefully--your choice will (most likely) affect your inferences. To see how, just look at the differences in the posterior summaries from this real data set using four different, but plausible, prior distributions. Note that Jeffreys' prior and the Cauchy(2.5) prior are two *defaults* recommended in the literature, yet the 90% credible interval for the latter is more than *twice as wide* as the former and the posterior mean is almost *twice as large*.

![An illustration that the choice of prior matters.](matters-ci.png)

And here's a working abstract:

> When facing small numbers of observations or ``rare events,'' political scientists often encounter separation, in which important explanatory variables perfectly predict binary events or non-events. In this situation, maximum likelihood provides implausible estimates and the researcher must incorporate some form of prior information into the model. The most sophisticated research uses Jeffreys' invariant prior to stabilize the estimates. While Jeffreys' prior has the advantage of being automatic, I show that it often provides too much prior information, producing point estimates that are too small and confidence intervals that are too narrow. To help researchers choose a more reasonable prior distribution, I introduce the concept of a partial prior predictive distribution and develop the tools required to compute the partial prior predictive distribution, estimate the subsequent model, and summarize the results.

If you have any comments or suggestions, please [open an issue](https://github.com/carlislerainey/priors-for-separatioin/issues) or just [e-mail](mailto:carlislerainey@gmail.com) me.