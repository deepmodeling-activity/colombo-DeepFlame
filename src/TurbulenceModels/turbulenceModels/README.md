# TurbulenceModels

## Using

- Add the following line into your 'controlDict':

  `libs ( "libdfTurbulenceModels.so" )`

- In 'constant/turbulenceProperties', e.g. for Sigma model, set:

  ```c++
  simulationType  LES;

  LES
  {
      LESModel Sigma;
      filter   simple;
  }
  ```
