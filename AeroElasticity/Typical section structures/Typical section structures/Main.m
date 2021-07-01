clear all
close all
clc

%% Get typical section parameters
Typical_section = GetTypicalSectionParameters;

%% Build structural matrices
[Ms,Ks] = GetTypicalSectionStructuralMatrices(Typical_section);