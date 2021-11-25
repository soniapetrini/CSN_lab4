# CSN_lab4: Non-Linear Regression on Dependency Trees
The goal of this lab is to fit __non-linear functions__ to _syntactic dependency trees_ from 10 different languages. 
Specifically, we focused on modeling the __degree second moment ⟨k2⟩__ distribution as a function of the __number of vertices (words)__ of a sentence, n. 
To do so, we implemented an ensemble of models, composed of a null model, some base models and their counterparts assuming an offset (+ models).
The Null model follows the distribution of uniformly random labelled trees, while the ensemble of base models includes a power-law, and exponential, and a logarithmic form.
