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

// get tree branch names

topology_branch_names = BranchName (bsrel_tree,-1);

matrices_defined          = 1;  
branch_length_expressions = {};
solve_these_for_lengths   = {};

// ----- apply nucleoitide biases

nucleotide_bias_settings ["apply"][""];
function apply (key, value) {
    ^key = value;
    return 0;
}

nonsyn_vals = bsrel_settings["nonsynvals"]
syn_vals = bsrel_settings["synvals"]

for (branch_id = 0; branch_id < Columns(topology_branch_names)-1; branch_id += 1) {
    branch_name = topology_branch_names[branch_id];
    full_name = "bsrel_tree.`branch_name`";

    PopulateModelMatrix("MGMatrix" + matrices_defined, nucCF, "syn", "", "nonsyn");
    ExecuteCommands("Model MGlocal = (MGMatrix" + matrices_defined + ", codon3x4, 0)");

    ExecuteCommands ("SetParameter (`full_name`, MODEL, MGlocal);");
    *(full_name + ".nonsyn") = nonsyn_vals[0];
    *(full_name + ".syn") = syn_vals[0];
}

codonCharacters = {{"A","C","G","T"}
			  			   {"3",GeneticCodeExclusions,"",""}};

SetDialogPrompt ("Save simulated data to:");
fprintf 		(PROMPT_FOR_FILE, CLEAR_FILE);
save_sims_to     = LAST_FILE_PATH;
IS_TREE_PRESENT_IN_DATA = 1;
DATAFILE_TREE = Format (bsrel_tree,1,1);

fprintf (stdout, "\n");


// only simulating a single codon (fourth parameter)
DataSet combinedSet = Simulate 	(bsrel_tree,codon3x4,codonCharacters,1,0);

// codons is defined in the settings file
for (it = 1; it < codons; it += 1) {
    for (branch_id = 0; branch_id < Columns(topology_branch_names)-1; branch_id += 1) {
        branch_name = topology_branch_names[branch_id];
        full_name = "bsrel_tree.`branch_name`";
        *(full_name + ".nonsyn") = nonsyn_vals[it];
        *(full_name + ".syn") = syn_vals[it];
    }
    DataSet sim = Simulate(bsrel_tree,codon3x4,codonCharacters,1,0);
    Dataset combinedSet = Concatenate(combinedSet, sim);
}

DataSetFilter simFilter = CreateFilter (combinedSet,1);
fName = save_sims_to + "." + it;
fprintf 		(fName, CLEAR_FILE,simFilter);
lfOut = save_sims_to + ".fit";
