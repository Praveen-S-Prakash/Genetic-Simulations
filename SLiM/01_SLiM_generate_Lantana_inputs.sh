#!/bin/bash

# SLiM - Lantana invasion simulation
#
# This script generates SLiM input files for 100 replicate simulations
# of an invasive population scenario.
#
# The simulation starts with an ancestral population of 100 individuals,
# allows the population to grow exponentially, and then splits it into
# three independently expanding invasive populations (p2, p3 and p4).
#
# The simulations use a 1 Mb nucleotide-based genome with a specified
# mutation rate and recombination landscape. Selfing rate, number of
# founding individuals, and VCF sample size are set as parameters below.
#
# Each replicate generates VCF files for the ancestral and three
# invasive populations, as well as a log containing population genetic
# summary statistics.

# Simulation parameters
run_on_cluster=1
Nrep=100
neme=_0_
Selfrate=0
Nsplit=10
VcfSampleN=50

pwdOUT=~/Slim_simulations_cluster/lantana_simulations/100_replicates/migration_10
scenar=scenario0


# Generate one SLiM input file for each replicate
for ((i=1; i<=Nrep; i++)); do

    # Initialize the genome
    #
    # A 1 Mb nucleotide-based genome is simulated with a Jukes-Cantor
    # mutation model and the specified recombination-rate landscape.

    echo -e "initialize() \n
    { \n
    initializeSLiMOptions(nucleotideBased=T); \n
    initializeAncestralNucleotides(randomNucleotides(1000000)); \n
    initializeMutationTypeNuc(\"m1\", 0.5, \"f\", 0.0); \n
    initializeGenomicElementType(\"g1\", m1, 1.0, mmJukesCantor(1e-7)); \n
    initializeGenomicElement(g1, 0, 999999); \n
    rates = c(1e-8, 0.5, 1e-8, 0.5, 1e-8, 0.5, 1e-8, 0.5, 1e-8, 0.5, 1e-8, 0.5, 1e-8, 0.5, 1e-8, 0.5, 1e-8, 0.5, 1e-8); \n
    ends = c(99999, 100000, 199999, 200000, 299999, 300000, 399999, 400000, 499999, 500000, 599999, 600000, 699999, 700000, 799999, 800000, 899999, 900000, 999999); \n
    initializeRecombinationRate(rates, ends); \n
    } \n" > $pwdOUT/$scenar/inputFile$neme$i.txt


    # Ancestral population and burn-in
    #
    # Start with 100 individuals and apply the specified selfing rate.
    # The population then grows exponentially until it reaches 10,000.

    echo -e "1 \nearly() \n
    { \n
    sim.addSubpop(\"p1\", 100); \n
    p1.setSelfingRate(Selfrate); \n
    } \n
    " >> $pwdOUT/$scenar/inputFile$neme$i.txt

    echo -e "1:20200 \nearly() \n
    { \n
    if (p1.individualCount < 10000) \n
    { \n
    newSizeP1 = asInteger(round(1.08^(sim.cycle - 0) * 100)); \n
    p1.setSubpopulationSize(newSizeP1); \n
    } \n
    }" >> $pwdOUT/$scenar/inputFile$neme$i.txt


    # Invasion: population splitting
    #
    # At generation 20,000, the ancestral population is split into
    # three independently founded invasive populations (p2, p3 and p4).
    # Each invasive population is founded by Nsplit individuals and
    # inherits the specified selfing rate.

    echo -e "20000 \nearly() \n
    { \n
    sim.addSubpopSplit(\"p2\", Nsplit, p1); \n
    p2.setSelfingRate(Selfrate); \n

    log = community.createLogFile(\"$pwdOUT/$scenar/sim_log_trial_SP${neme}P10000_M10_S50_$i.txt\", logInterval=10); \n
    log.addCycle(); \n
    log.addCustomColumn(\"FST_p1_p2\", \"calcFST(p1.genomes, p2.genomes);\"); \n
    log.addCustomColumn(\"FST_p1_p3\", \"calcFST(p1.genomes, p3.genomes);\"); \n
    log.addCustomColumn(\"FST_p1_p4\", \"calcFST(p1.genomes, p4.genomes);\"); \n
    log.addCustomColumn(\"FST_p2_p3\", \"calcFST(p2.genomes, p3.genomes);\"); \n
    log.addCustomColumn(\"FST_p2_p4\", \"calcFST(p2.genomes, p4.genomes);\"); \n
    log.addCustomColumn(\"FST_p3_p4\", \"calcFST(p3.genomes, p4.genomes);\"); \n
    log.addCustomColumn(\"Mean_heterozygosity_p1\", \"calcHeterozygosity(p1.genomes);\"); \n
    log.addCustomColumn(\"Mean_heterozygosity_p2\", \"calcHeterozygosity(p2.genomes);\"); \n
    log.addCustomColumn(\"Mean_heterozygosity_p3\", \"calcHeterozygosity(p3.genomes);\"); \n
    log.addCustomColumn(\"Mean_heterozygosity_p4\", \"calcHeterozygosity(p4.genomes);\"); \n
    log.addCustomColumn(\"WattersonsTheta_p1\", \"calcWattersonsTheta(p1.genomes);\"); \n
    log.addCustomColumn(\"WattersonsTheta_p2\", \"calcWattersonsTheta(p2.genomes);\"); \n
    log.addCustomColumn(\"WattersonsTheta_p3\", \"calcWattersonsTheta(p3.genomes);\"); \n
    log.addCustomColumn(\"WattersonsTheta_p4\", \"calcWattersonsTheta(p4.genomes);\"); \n
    }" >> $pwdOUT/$scenar/inputFile$neme$i.txt

    echo -e "20000 \nearly() \n
    { \n
    sim.addSubpopSplit(\"p3\", Nsplit, p1); \n
    p3.setSelfingRate(Selfrate); \n
    } \n
    " >> $pwdOUT/$scenar/inputFile$neme$i.txt

    echo -e "20000 \nearly() \n
    { \n
    sim.addSubpopSplit(\"p4\", Nsplit, p1); \n
    p4.setSelfingRate(Selfrate); \n
    } \n
    " >> $pwdOUT/$scenar/inputFile$neme$i.txt


    # Invasive population growth
    #
    # After founding, p2, p3 and p4 independently expand from
    # the founding population toward a maximum size of 10,000.

    echo -e "20000:20200 \nearly() \n
    { \n
    if (p2.individualCount < 10000) \n
    { \n
    newSizeP2 = asInteger(round(1.08^(sim.cycle - 19999) * 10)); \n
    p2.setSubpopulationSize(newSizeP2); \n
    } \n

    if (p3.individualCount < 10000) \n
    { \n
    newSizeP3 = asInteger(round(1.08^(sim.cycle - 19999) * 10)); \n
    p3.setSubpopulationSize(newSizeP3); \n
    } \n

    if (p4.individualCount < 10000) \n
    { \n
    newSizeP4 = asInteger(round(1.08^(sim.cycle - 19999) * 10)); \n
    p4.setSubpopulationSize(newSizeP4); \n
    } \n
    }" >> $pwdOUT/$scenar/inputFile$neme$i.txt


    # Sample genomes and save VCF files
    #
    # At the end of the simulation, sample genomes from each population
    # and export them as VCF files for downstream analyses.

    echo -e "20200 late() \n
    { \n
    p1.outputVCFSample(VcfSampleN, filePath=(\"$pwdOUT/$scenar/p1_SP${neme}P10000_M10_S50_$i.vcf\")); \n
    p2.outputVCFSample(VcfSampleN, filePath=(\"$pwdOUT/$scenar/p2_SP${neme}P10000_M10_S50_$i.vcf\")); \n
    p3.outputVCFSample(VcfSampleN, filePath=(\"$pwdOUT/$scenar/p3_SP${neme}P10000_M10_S50_$i.vcf\")); \n
    p4.outputVCFSample(VcfSampleN, filePath=(\"$pwdOUT/$scenar/p4_SP${neme}P10000_M10_S0_$i.vcf\")); \n
    }" >> $pwdOUT/$scenar/inputFile$neme$i.txt

done
