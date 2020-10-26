M = read_microenvironment( '/Users/dom/PhysiCell/haskell/PhysiCell/output/initial_microenvironment0.mat' );

titles{1} = 'cells';
titles{2} = 'blood vessels';
titles{3} = 'oxygen';
plot_microenvironment( M , titles );