# glmmr
This package provides a suite of design analysis tools for any study design that can be represented or modelled as a generalised linear mixed model. These studies include cluster randomised trials, cohort studies, spatial and spatio-temporal epidemiological models, and split-plot designs. The aim of the package is to provide flexible access to various methods like power calculation, either by simulaton or approximation, or identification of optimal designs. 

We use the R6 class system to represent four types of object: mean function, covariance, a design, and a design space. Each of these objects takes different inputs and produces and stores different components of the GLMM model. The intention of representing models in this modular fashion is to permit users to change individual components, such as a single parameter or covariance specification, and have all linked components update seamlessly. 

A complete vignetted will be added here later.
