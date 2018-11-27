# RBF-network
Design an RBF network without using any existing machine learning library, including libraries for the k-means algorithm. 

1. The goal is to design an RBF network g(x) with 20 centers. Run the k-means
algorithm for 10 centers for class C1. Set these as centers c1, . . . , c10, and sketch them. Run the k-means
algorithm for 10 centers for class C−1. Set these as centers c11, . . . , c20. Sketch these as well, but use
different markers compared to ones you used for centers of class C1.
2. Now, run the perceptron training algorithm to determine the weights ω1, . . . , ω20 and the bias θ. If you
are doing everything correctly, your PTA should converge, and you should be able to separate the two
classes perfectly. Provide a rough sketch of the corresponding decision boundary {x : g(x) = 0}.
3.  Repeat 1. and 2. for the case of a total of 4 centers. Again use half of the centers for one class, and the
other half for the other. Comment on the differences (if any).
