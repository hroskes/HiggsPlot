vector<TString> branches;
vector<vector<void*> > variables;

int findIndex(TString branch);
unsigned int getMaxSize();

int SetBranchAddress(TTree *tree,TString branch,void *variable)
{
    int result = tree->SetBranchAddress(branch,variable);
    if (result < 0)
        return result;
    if (findIndex(branch) == -1)
    {
        branches.push_back(branch);
        vector<void*> v;
        variables.push_back(v);
    }
    variables[findIndex(branch)].push_back(variable);
    return result;
}

void GetEntry(TTree *tree,Long64_t entry = 0,Int_t getall = 0)
{
    for (unsigned int i = 0; i < getMaxSize(); i++)
    {
        for (unsigned int j = 0; j < branches.size(); j++)
        {
            if (variables[j].size() == 1) continue;
            tree->SetBranchAddress(branches[j],variables[j][i]);
        }
        tree->GetEntry(entry,getall);
    }
}

int findIndex(TString branch)
{
    for (unsigned int i = 0; i < branches.size(); i++)
        if (branches[i] == branch)
            return i;
    return -1;
}

unsigned int getMaxSize()
{
    unsigned int result = 0;
    for (unsigned int i = 0; i < branches.size(); i++)
        if (variables[i].size() > result)
            result = variables[i].size();
    return result;
}

//Try to set it to branches from the array, in order, until one works.
//Separate functions for 2 or 3, for convenience

int SetBranchAddress(TTree *tree,int nBranches,TString *trybranches,void *variable)
{
    int result = -6;
    for (int i = 0; i < nBranches; i++)
    {
        result = SetBranchAddress(tree,trybranches[i],variable);
        if (result >= 0) break;
    }
    return result;
}

int SetBranchAddress(TTree *tree,TString branch1,TString branch2,void *variable)
{
    TString trybranches[2] = {branch1,branch2};
    return SetBranchAddress(tree,2,trybranches,variable);
}

int SetBranchAddress(TTree *tree,TString branch1,TString branch2,TString branch3,void *variable)
{
    TString trybranches[3] = {branch1,branch2,branch3};
    return SetBranchAddress(tree,3,trybranches,variable);
}
