README
======

This is a Github repository with a set of Matlab scripts for reverse shooting. They originally come from Chad Jones' paper [Life and Growth](https://web.stanford.edu/~chadj/LifeandGrowthJPE2016.pdf), and were then modified for use on [Existential Risk and Growth](https://leopoldaschenbrenner.github.io/xriskandgrowth/ExistentialRiskAndGrowth050.pdf). I then cleaned it up and added comments, starting with some notes by Ben Snodin. 

The code contains two parts, one which carries out some vanilla reverse shooting, and another which adds shocks and accelerations. I have only cleaned out the first part, and don't particularly intend to engage with the second.

## Vanilla reverse shooting
Files: Calibrate.m, getsteadystate.m, transit1dx.m, solvetransition.m, getells0.m, Transition.m 

How to use: Change the "matlabminiscriptspath" variable on Transition.m, and then run Transition.m

Modify: Change the hard-coded values and equations in the vanilla reverse shooting files to your desired equations and values. Change the "matlabminiscriptspath" variable on Calibration.m as well.

## Shocks and accelerations
Files: Acceleration.m, TransitoryShock.m, FindAlternativePath.m 

How to use and modify: Here be dragons. Don't. Leave that to less innocent souls. 

## What each file does

### getsteadystate.m
Give the long-run values for important variables for some important cases (proposition 2 and proposition 5)

### transit1dx.m
A function representing a system of 6 ODEs for the time evolution of the model. Specifically, I think it gives the changes in the time-dependent variables corresponding to the current values and step-size

### solvetransition.m
A function that computes the path for the system over the specified time interval. It finds the "best" s_0 and l_0 given the other parameters provided, and returns the transition corresponding to that. It uses getells0.m to get the best s_0 and l_0

### getells0.m
Returns values of s_o and l_o such that d(l)_0 - g_s* and d(s)_0 - g_s* and d(sigma)_0 are minimized, given all other model parameters (note that you don't need a transition path for this).  d(some variable) is that variable with a ^ on top, in the paper. Represents the derivative

### Calibrate.m
prints values for various values that represent a calibration that's specified in Appendix B1 of the paper. The objective function is based on the transition path. Note that, for each proposed choice of these variables, s_0 and l_0 are found using getells0.m, which is called by solvetransition.m as it gets the transition path.

### Transition.m 
Finds the transition path using solvetransition.m using the values from Calibrate.m for the various variables (note that these values are hard-coded) Plots some variables for the path (I think I recognise them from the paper)

### Acceleration.m
Finds a transition path using solvetransition.m which represents the baseline case. Then uses FindAlternativePath to find the path with accelerated growth. FindAlternativePath.m works by generating paths with accelerated growth with various d_0 and N_0, choosing one that at some point passes through the same state as teh baseline case does at some time. Plots some variables from both paths (I think I recognise them from the patper)

### TransitoryShock.m
As Acceleration.m, but for the transitory shock rather than the acceleration case. 

### FindAlternativePath.m 
Returns values of delta_0 and N_0 such that the associated path minimizes the following objective function "find the time t_0 on the path where delta_(t_0) is close to a given target. Then get the weighted sum of squared deviations between the state variables and the target state variables. 

## General notes and gotchas
Matlab doesn't seem to have a return keyword; instead return variables are specified with the function definition.

I.e., what would be 

```
someFunctionName = function(some variables){
  do something here
  return(results)
}
```

in any reasonable programming language, is instead in matlab

```
[result1, result2, ...] = someFunctionName(some variables)
  do something here
end
```

In functions which take a whole file (and which share the name with that file) it isn't necessary to write the end keyword, but I have because it's less confusing.

This code may not work on Octave, an open source matlab clone, and didn't as of 4/Sept/2020, because optimoptions and optimset have [not been exported](https://wiki.octave.org/Optimization_package) to Octave yet.

In matlab scripts (i.e., files with a .m extension, rather than inputs into the matlab terminal), function definitions must be placed at the end of a file, even if they are used before.

Logs are created on the directory returned by the pwd (print working directory) matlab function. 
