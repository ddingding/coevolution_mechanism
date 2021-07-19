# config file

# antitoxin single mutant library in presence of wild-type toxin is labelled as 4, and replicates are labelled as 1
# or 2 following the dot, the dash at the end indicates the timepoint the samplew as taken.
expToInd = {
    "4.1_0"  : "TAAGGCGA",
    "4.2_0"  : "CGTACTAG",
    "4.1_600": "GGACTCCT",
    "4.2_600": "TAGGCATG",
}

expToTemplate = {"4.1_0": "pard", "4.2_0": "pard", "4.1_600": "pard", "4.2_600": "pard"}

exp_to_od = {
    "4.1_0"  : 0.094,
    "4.2_0"  : 0.087,
    "4.1_600": 86.3927395,
    "4.2_600": 59.95732414,
}

tc_to_cond = {"41": "AT singles", "42": "AT singles"}
