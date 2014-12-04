RequireVersion ("2.1320130313");
VERBOSITY_LEVEL				= 0;
skipCodeSelectionStep 		= 0;

SetDialogPrompt      ("Simulation settings:");
ExecuteAFile         (PROMPT_FOR_FILE);

LoadFunctionLibrary("chooseGeneticCode", {"0": "Universal"});
LoadFunctionLibrary("BranchSiteTemplate");
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

//global nonsyn = nonsynvals[0][0];
//global syn = synvals[0][0];

for (branch_id = 0; branch_id < Columns(topology_branch_names)-1; branch_id += 1) {
    branch_name = topology_branch_names[branch_id];
    full_name = "bsrel_tree.`branch_name`";

    PopulateModelMatrix("MGMatrix" + matrices_defined, nucCF, "syn", "", "nonsyn");
    ExecuteCommands("Model MGlocal = (MGMatrix" + matrices_defined + ", codon3x4, 0)");

    ExecuteCommands ("SetParameter (`full_name`, MODEL, MGlocal);");

    length_info = (bsrel_settings[branch_name])["length"];
    length = length_info[0];
    *(full_name + ".nonsyn") = nonsynvals[0][0] * length;
    *(full_name + ".syn") = synvals[0][0] * length;
    //fprintf(stdout, "" + length[0] + "\n");
    //fprintf(stdout, "" + nonsynvals[0][0] * length + "\n");
    //fprintf(stdout, "" + synvals[0][0] * length + "\n");
    //*(full_name + ".nonsyn") = nonsynvals[0][0];
    //*(full_name + ".syn") = synvals[0][0];
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
        //*(full_name + ".nonsyn") = nonsynvals[0][it];
        //*(full_name + ".syn") = synvals[0][it];
        //nonsyn = nonsynvals[0][it];
        //syn = synvals[0][it];
        length_info = (bsrel_settings[branch_name])["length"];
        length = length_info[0];
        *(full_name + ".nonsyn") = nonsynvals[0][0] * length;
        *(full_name + ".syn") = synvals[0][0] * length;
    }
    DataSet sim = Simulate(bsrel_tree,codon3x4,codonCharacters,1,0);
    DataSet combinedSet = Concatenate(combinedSet, sim);

    /*
    lfOut = save_sims_to + "." + it + ".fit";
    DataSetFilter newFilter = CreateFilter(sim, 3, "", "", GeneticCodeExclusions);
    LikelihoodFunction sim_LF = (newFilter, bsrel_tree);
    LIKELIHOOD_FUNCTION_OUTPUT = 7;
    fprintf(lfOut, CLEAR_FILE, sim_LF);
    LIKELIHOOD_FUNCTION_OUTPUT = 2;
    */
}

DataSetFilter simFilter = CreateFilter (combinedSet,1);
fName = save_sims_to;
fprintf 		(fName, CLEAR_FILE,simFilter);
