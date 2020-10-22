# How to write a .json config file for the simulation

## Quick summary of the JSON format
JSON files are structures with key-value pairs. Keys can only be strings, values can be strings, numbers (real-valued), lists (of objects of the same type) or other JSON objects.

Comments are supported by some JSON specifications but are not part of the original JSON format.  
**Comments are *not* supported by the C++ JSON parser (nlohmann-json) this project is using**.

```jsonc
{
    "key-string": "value",
    "key-number": 42,
    "key-list": [1,1,2,3,5,8,13],
    "key-object": {
        "key-string": "value",
        "key-number": 42,
        "key-list": [1,1,2,3,5,8,13],
        "key-object": {}
    }
}
```

---

## Structure of the configuration file
The config file for a simulation case is structured into four main parts. All values are in SI units unless otherwise stated.

```jsonc
{
    // base config goes here,
    "starting_conditions": {
        // ...
        "cloud" : {
            // ...
        }
    },
    "potentials": [
        {
            // ...
        },
        // ...
    ],
    "timepoints": [
        {
            // ...
        },
        // ...
    ]
}
```

---

### Base config
> | key | explanation | expected datatype | optional |
> | --- | ----------- | ----------------- | -------- |
> | `runtime` | the maximum duration | numeric | :no_entry_sign: |
> | `mass` | the mass of the atomic species | numeric | :no_entry_sign: |
> | `cell_structure` | the number of DSMC cells in x,y,z direction | list with three integers | :no_entry_sign: |

---

### Initial conditions
> | key | explanation | expected datatype | optional |
> | --- | ----------- | ----------------- | -------- |
> | `atom_number` | the number of atoms to start with | integer | :no_entry_sign: |
> | `temperature` | the initial temperature of the cloud | numeric | :no_entry_sign: |
> | `cloud` | the initial shape of the cloud | json object | :no_entry_sign: |
> | `scattering_length` | the s-wave scattering length of the atoms | numeric | :no_entry_sign:
> | `losses` | the background and 2-/3-body loss rates, defaults all to 0 if omitted | json object | :white_check_mark: |

#### Cloud
> | key | explanation | expected datatype | optional |
> | --- | ----------- | ----------------- | -------- |
> | `shape` | cloud density distribution | string: either `"uniform"` or `"harmonic"` | :no_entry_sign: |
> | `side_length` | if shape is 'uniform' -> side length of the cloud | numeric | :exclamation: mutually exclusive with `radius` | 
> | `radius` | if shape is 'harmonic' -> RMS radius of the cloud | numeric | :exclamation: mutually exclusive with `side_length`
> | `velocity_distribution` | max-well boltzmann distribution or same speed for all | string: either `"maxwell-boltzmann"` (default if omitted) or `"delta"` | :white_check_mark: |

#### Losses
Possible entries are `K1` for background losses (in 1/s), `K2` for two-body inelastic losses (in m^3/s) and `K3` for three-body inelastic losses (in m^6/s).  
All are optional and will default to 0 if omitted.

---

### Potentials
Each entry in this list contains:

> | key | explanation | expected datatype | optional |
> | --- | ----------- | ----------------- | -------- |
> | `type` | the type of potential, as registered in `main.cpp` | string | :no_entry_sign: |
> | `id` | a unique (but arbitrary) identifier for the potential | string | :no_entry_sign: |
> | `active` | if the potential is active from the beginning, defaults to `true` | boolean | :white_check_mark: |
> | `parameters` | dictionary of parameters for this potential | json object | :grey_question: depends on `type` |

E.g. the gravity potential doesn't have any parameters, therefore the `parameters` entry can be omitted in that case.  
Otherwise the parameters are named in the definition for the respective potential. Parameters (if existing) **must** have numeric values.

---

### Timepoints
:exclamation: This whole entry is **optional** :exclamation:  
Specifies the time evolution of the potential parameters (and active setting) as well as for the scattering length and the loss parameters.  
Each entry in this list **must** contain a `time` parameter which specifies the end of the time-period (absolute), meaning entries with 1.5, 3.5 specify the times 1.5s and 3.5s, **not** intervals of 1.5s and 3.5s, so in this case the intervals would be 1.5 and 2 seconds long.  

> | key | explanation | expected datatype | optional |
> | --- | ----------- | ----------------- | -------- |
> | `time` | see above | numeric | :no_entry_sign: |
> | `potentials` | a list of changes for potential parameters | list of json objects (see below) | :white_check_mark: |
> | `scattering_length` | change of the scattering length | json object | :white_check_mark: |
> | `losses` | changes to the loss parameters `K{1-3}` | json objects | :white_check_mark: |

Each change is defined by a json object:

> | key | explanation | expected datatype | optional |
> | --- | ----------- | ----------------- | -------- |
> | `start` | initial value, defaults to the starting value at t=0 or the final value at the last timepoint | numeric | :white_check_mark: |
> | `end` | final value at the end of the interval | numeric | :no_entry_sign: |
> | `model` | functional form of the change, as named in the `Models` factory singleton (`model-factory.hpp`) | string | :no_entry_sign: |
> | `parameters` | additional parameters for the function (e.g. time constant for an exponential) | json object | :grey_question: depends on `model` |

To change a potential parameter, the potential needs to be identified by its `id`.  
It is also possible to switch a potential from active to inactive or vice versa.  
A list entry of the `potentials` parameter then has the following entries:

> | key | explanation | expected datatype | optional |
> | --- | ----------- | ----------------- | -------- |
> | `id` | the unique identifier given to the potential before | string | :no_entry_sign: |
> | `active` | if the potential should be active/inactive **after** the specified time | boolean | :white_check_mark:
> | parameter name | the change for this specific parameter | json object for a change (see above) | :white_check_mark: |

---

# Full example of a configuration file with comments

```jsonc
{
    "runtime": 10, // run for 10s
    "mass": 42, // each atoms weighs 42kg
    "cell_structure": [7,7,7], // 7 DSMC cells in each direction (343 in total)
    "starting_conditions": {
        "atom_number": 1000000000, // start with 1 billion atoms
        "temperature": 1e-6, // 1uK
        "cloud": {
            "shape": "uniform", // initialize atoms in a cube
            "side_length": 1 // with side length 1m
        },
        "scattering_length": 1e-3, // a = 1mm
        "losses": {
            "K1": 0.016667, // lifetime of 1 minute
            "K3": 1e-40 // three body loss rate 1e-40 m^6/s
            // K2 will be 0 by default in this case
        }
    },
    "potentials": [
        {   // first potential
            "type": "some_potential", // potential named "some_potential" in main.cpp
            "id": "funny id for this potential", // please don't actually name your potentials like this
            "parameters": {
                "some_parameter": 3,
                "other_parameter": 5
            }
            // active defaults to true
        },
        {   // second potential
            "type": "gravity", // I'd recommend using unmistakeable names
            "id": "gravity is not a force",
            // no parameters necessary
            "active": false // gravity will not have any effect
        }
    ],
    "timepoints": [
        {
            "time": 1, // first interval ends at t = 1s
            "losses": {
                "K2": {
                    "end": 1e-20, // after 1s, the two-body loss rate will be 1e-20 m^3/s
                    "model": "switch" // the change will be instantaneous, at t = 0.99, K2 = 0, at t = 1, K2 = 1e-20
                },
                "K3": {
                    "end": 2e-40,
                    "model": "exp",
                    "parameters": {
                        "tc": 2 // time constant
                    }
                }
            },
            "potentials": [
                {
                    "id": "gravity is not a force",
                    "active": true // gravity will turn on after 1s
                }
            ]
        },
        {
            "time": 4, // second interval ends at t = 4s
            "potentials": [
                {
                    "id": "funny id for this potential",
                    "some_parameter": {
                        "end": 7,
                        "model": "lin"
                    } // "some_parameter" will change linearly from 3 at t=0 to 7 at t=4
                },
                {
                    "id": "gravity is not a force",
                    "active": false // gravity will turn back off
                }
            ],
            "losses": {
                "K3": {
                    "start": 1e-40,
                    "end": 4e-40,
                    "model": "lin"
                } // K3 will change instantaneously from 2e-40 back to 1e-40 at t=1s and then evolves linearly to 4e-40 at t=4s
            }
        },
        {
            "time": 8.5, // third interval ends at t = 8.5s
            "potentials": [
                {
                    "id": "funny id for this potential",
                    "other_parameter": {
                        "end": 0,
                        "model": "exp",
                        "parameters": {
                            "tc": 0.1
                        }
                    }, // "other_parameter" will decay exponentially from 5 at t=0 to 0 at t=8.5 with a short time constant of 0.1s
                    "active": false // "funny id for this potential" will turn off completely after 8.5 seconds
                },
                {
                    "id": "gravity is not a force",
                    "active": true // gravity will turn back on
                }
            ]
        }
    ]
}
```

The different parameters evolve completely independent from each other.  
There can be arbitrarily many timepoints.
