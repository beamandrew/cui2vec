# cui2vec

This repo contains the code associated with the following paper (under review):

> Beam, A.L., Kompa, B., Schmaltz, A., Fried, I., Griffin, W, Palmer, N.P., Shi, X., Cai, T., and Kohane, I.S.,, 2019. Clinical Concept Embeddings Learned from Massive Sources of Multimodal Medical Data. arXiv preprint arXiv:1804.01486.

# Overview

This repo contains the R pacakge `cui2vec`, which provides code for fitting embeddings to your own co-occurrence data in the manner presented in the above paper. The package can be installed locally from source. An overview of usage is provided in the following HTML vignette, which can be viewed in your browser:

[vignettes/rendered/2019_07_31/cui2vecWorkflow.html](vignettes/rendered/2019_07_31/cui2vecWorkflow.html).

Additional information on each of the public functions can be accessed in the standard way (e.g., ```?cui2vec::construct_word2vec_embedding```).

Data agreements prevent us from releasing all of our original source data, but upon acceptance, we will release our embeddings at the following URL: TBD.
