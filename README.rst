.. role:: math(raw)
   :format: html latex
..

**Gaussian-Process-Model-based**

**System-Identification**

**Toolbox for Matlab**

**Version 1.2.2**

Martin Stepančič and Juš Kocijan

.

Introduction
============

The idea of this toolbox is to facilitate dynamic systems identification
with Gaussian-process (GP) models. The presented toolbox is continuously
developing and is put together with hope to be useful as a springboard
for the modelling of dynamic systems with GP models.

The GP model belongs to the class of black-box models. GP modelling
differs from most other black-box identification approaches in that it
does not try to approximate the modelled system by fitting the
parameters of the selected basis functions, but rather it searches for
the relationship among the measured data. The model is composed of
input-output data that describes the behaviour of the modelled system
and the covariance function that describes the relation with respect to
the input-output data. The prediction of the GP model output is given as
a normal distribution, expressed in terms of the mean and the variance.
The mean value represents the most likely output, and the variance can
be interpreted as a measure of its confidence.

System identification is composed of methods to build mathematical
models of dynamic systems from measured data. It is one of the
scientific pillars used for dynamic-systems analysis and control design.
The identification of a dynamic system means that we are looking for a
relationship between past observations and future outputs.
Identification can be interpreted as the concatenation of a mapping from
measured data to a regression vector, followed by a nonlinear mapping
from the regression vector to the output space. Various machine-learning
methods and statistical methods are employed to determine the nonlinear
mapping from the regression vector to the output space. One of the
possible methods for a description of the nonlinear mapping used in
identification is GP models. It is straightforward to employ GP models
for the discrete-time modelling of dynamic systems within the
prediction-error framework.

Many dynamic systems are often considered as complex; however,
simplified input-output behaviour representations are sufficient for
certain purposes, e.g., feedback control design, prediction models for
supervisory control, etc.

| More on the topic of system identification with GP models and the use
  of this models for control design can be found in the book:
| Juš Kocijan (2016) Modelling and Control of Dynamic Systems Using
  Gaussian Process Models, Springer.

GP-Model-based System-Identification Toolbox for Matlab
=======================================================

Prerequisites
-------------

As this toolbox is intended to use within Matlab environment the user
should have Matlab installed. It works on Matlab 7 and later, but there
should be no problems using the toolbox on previous versions of Matlab,
e.g., 6 or 5.

It is also assumed that the GPML toolbox [1]_, general purpose GP
modelling toolbox for Matlab, is installed. The GP-model-based
system-identification toolbox serves as upgrade to GPML toolbox.

The user should posses some familiarity with the Matlab structure and
programming.

Installing GPdyn toolbox
------------------------

Unzip the file GPdyn into chosen directory and add path, with
subdirectories, to Matlab path.

Overview of the GPdyn toolbox
-----------------------------

GPdyn files are contained in several directories, depending on their
purpose:

training functions,
    used for training GP models of dynamic systems;

GP-model evaluation functions,
    used for simulating the dynamic GP model;

LMGP-model evaluation functions,
    which are used when modelling and simulating the system with a GP
    model with incorporated local models (LMGP model);

utilities functions,
    that are various support functions;

demo functions,
    which demonstrate the use of the toolbox for identification of
    dynamic systems.

The list of included functions, demos and one model is given in
following tables.

| \|l\|l\|
|  & GP-model training of ARX model
|  & GP-model training of OE model
|  & - finding initial values of hyperparameters with random search
|  & minimize a multivariate function using differential evolution
|  & GP model simulation without the propagation of uncertainty
|  & GP model simulation with Monte Carlo approximation
|  & GP model simulation with analytical approximation of statistical
| & moments with a Taylor expansion for the squared exponential
| & covariance function
|  & GP model simulation with exact matching of statistical moments
| &for the squared exponential covariance function
|  & GP model simulation with exact matching of statistical moments
| &for the linear covariance function
|  & multi-step-ahead prediction of GP model without
| &the propagation of uncertainty
|  & modified version of GP rutine from the GPML toolbox
|  & creates samples of mixture components
|  & GP model prediction with stochastic inputs for
| &the squared exponential covariance function with Taylor expansion
|  & GP model prediction with stochastic inputs for
| &the linear covariance function
|  & GP model prediction with stochastic inputs for
| &the squared exponential covariance function

\|l\|l\|c\|

| 
|  & LMGP model simulation without the propagation of uncertainty
|  & LMGP model simulation with Monte Carlo approximation
|  & LMGP model training
|  & - LMGP model prediction
| & - data likelihood and its derivatives

| 
|  & adding white noise to noise-free simulation results
|  & construction of the input regressors
| & from system’s input signals
|  & method to evaluate covariance, mean and likelihood functions
|  & calculates negative log marginal likelihood
|  & the method for the lag-space selection, based on Lipschitz
  quotients
|  & checking of the parameters match
|  & performance measures
|  & testing sampled probability distributions
|  & plot results (output and error) of the GP model prediction
|  & plot error of the GP model prediction
|  & plot output of the GP model prediction
|  & preprocessing of data
|  & postprocessing of data
|  & postprocessing of predicted variance
|  & generating pseudo-random binary signal
|  & generating pseudo-random signal

\|l\|l\|

| 

| 

| & present the system used in demos
|  & generate data for the identification and validation
| & of the GP model
|  & normalization of input and output data
|  & training of the GP model
|  & validation with simulation of the GP model
|  & generate data for the identification and validation
| & of the LMGP model
|  & training of the LMGP model
|  & simulation of the LMGP model
|  & system simulation
|  & obtaining system’s derivatives
|  & identification of system’s local models

How to use this toolbox
-----------------------

Demos
~~~~~

| A simple nonlinear dynamic system is used to demonstrate the
  identification and simulation of the GP models:

  .. math:: y(k+1) = \frac{y(k)}{1+y^2(k)} + u^3(k) \label{eq:narendra}

  The system was used as an example of dynamic system identification
  with artificial neural networks in:
| K.S. Narendra and K. Parthasarathy. Identification and Control of
  Dynamical Systems Using Neural Networks, IEEE Transactions on Neural
  Networks, Vol.1 No. 1, 4–27, 1990.

demo\_example\_present,
    presents this system.

Following three demos present the identification of dynamic systems with
the GP model:

demo\_example\_gp\_data,
    which presents how to obtain and assemble data for identification;

demo\_example\_gp\_norm,
    which shows how to normalise input and output data for training;

demo\_example\_gp\_training,
    which demonstrates the identification with a GP model;

demo\_example\_gp\_simulation,
    which shows how to simulate the GP model.

The use of the GP model with incorporated local models is presented with
demos:

demo\_example\_lmgp\_data,
    which presents how to obtain and assemble data for identification;

demo\_example\_lmgp\_training,
    which demonstrates the training (=identifying) the LMGP model;

demo\_example\_lmgp\_simulation,
    which shows how to simulate the LMGP model.

Acknowledgements
~~~~~~~~~~~~~~~~

We would like to thank all past, present and future contributors to this
toolbox.

.. [1]
   It can be obtained from *http://www.gaussianprocess.org/gpml*.
