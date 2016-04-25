clear all; close all; clc;

load('constrOrig')
load('constrOrig_trunc4')
load('constrOrig_trunc3')

load('constrZim')
load('constrZim_trunc4')
load('constrZim_trunc3')


normZ_C=norm(C)
normO=norm(Corig)
normZ4=norm(Czim_trunc4)
normO4=norm(Corig_trunc4)
normZ3=norm(Czim_trunc3)
normO3=norm(Corig_trunc3)


