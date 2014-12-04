RequireVersion ("2.1320130313");
VERBOSITY_LEVEL				= 0;
skipCodeSelectionStep 		= 0;

SetDialogPrompt      ("Simulation settings:");
ExecuteAFile         (PROMPT_FOR_FILE);

LoadFunctionLibrary("chooseGeneticCode");
LoadFunctionLibrary("GrabBag");
LoadFunctionLibrary("dSdNTreeTools");
LoadFunctionLibrary("CF3x4");
LoadFunctionLibrary("TreeTools");

nucCF						= CF3x4	(nuc3x4, GeneticCodeExclusions);
codon3x4					= BuildCodonFrequencies (nucCF);

// Generate Distributions


// Generate Models

// Simulate sequences one site at a time, one sequence at a time, maintaining
// the input model pattern order
sequences = {};
for (seq = 
