Guan and Haran (2018)

## 9/2 Notes

### Notes

- The method seems to perform much much much faster than others
- Pretty straigtforward, still not sure about the estimation step and what that looks like in reality. It looks like an implementation would just be a lot of CLA style stuff which I can do. 
- Not really sure how to add a other mixed effect to the model and what considerations would be required. Seems like I could just toss it in and add it to the column space that the spatial mixed effects need to be orthogonal to. Or maybe I should also ensure the new mixed effects and the fixed effects are orthogonal too? Ask Alicia about that. 
- Line from paper: Developing extension of this methodology to spatial-temporal and multivariate spatial processes may provide fruitful avenues for future research. (hopefully I wouldn't be stepping on their toes lol.)

### Questions

- What does it mean to adjust a model after the fact to fix the coverage? This paper just says the adjusted coefficients fit well while the normal coeffiecients have low coverage (40%). I belive coverage is how well the posterior predictive intervals match up.
    - They cite this paper as an explanation of how to do this adjustument: Restricted spatial regression in practice: Geostatistical models, confounding, and robustness under model misspecication. (Hanks et al 2015)
- Not fully sure what the diff between the FRP and RRP and RSR are. 
- How is this meaningfully differnet from the Hughes and Haran (2013) method? I think its just that this can be used in point processes because it doesn't rely on the Moran prior, which doesn't work well for point processes.
- Why do they use this fancy random projections using the Nystrom method. How is this different from a normal PCA/SVD? Maybe its just a quick way of getting the first couple eigenvalues and vectors for a positive definite matrix?

