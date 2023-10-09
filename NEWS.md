# changes in version 1.2.1
* Enhancements
  - phase2.TTE() replaces and combines the functionality of Optimal.KJ() and 
  Optimal.rKJ() in previous versions. In addition to the optimal method, this 
  new function also provides designs for single-arm phase II survival trials using 
  the minmax and admissible methods.
  
  - Sim() replaces and combines the functionality of Sim_KJ() and Sim_rKJ() in 
  previous versions.
  
  - Examples about how to use phase2.TTE() and Sim() are provided in the 
  vignette and help documentations.


# changes in version 1.1.5

* Enhancements
  - for Optimal.KJ() and Optimal.rKJ(), now users can choose to specify a 
  minimum value for the early stopping probability under H0 using the argument 
  prStop. Examples about how to implement this new feature are provided in the 
  vignettes and help documentations. 
