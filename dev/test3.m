clear; matlabrc; clc; close all; rng(1);
import = @(x) addpath(genpath(x));
import('../data')

num_lmks = 100;
bennu = Bennu(num_lmks);

bennu.drawBody()
rotate3d on

bennu.drawLmks(true,'MarkerSize',20)