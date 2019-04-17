# Undergraduate Thesis Code

The thesis title is: "Study on the observability of the process H->J/ψγ->μμ with CMS experiment in HL-LHC".

## Repository structure
```bash
.
├── bkg
│   ├── pythia8gen
│   └── pythia8gensim
├── examples
│   ├── pythia8gen
│   └── pythia8gensim
├── multicpu
│   ├── bkg
│   │   ├── pythia8gen
│   │   └── pythia8gensim
│   ├── examples
│   │   ├── pythia8gen
│   │   └── pythia8gensim
│   ├── generator*
│   ├── generator*
│   ├── generator*
│   ├── generator*
│   └── signal
└── signal
    ├── pythia8gen
    └── pythia8gensim
```

```bkg``` directory contains scripts for production of J/ψ events using Pythia8 and simulation of detector response using Delphes.
```signal``` directory contains similar scripts dedicated to the process H->J/ψγ->μμ.
```multicpu``` directory contains parallelized implementations for both signal and background.
