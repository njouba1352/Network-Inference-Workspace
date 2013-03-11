import os, sys
from datetime import datetime
from DataStore import *
from ReadData import *
from JobManager import *
from Network import *
from Generate_Grid import *
from ReadConfig import *
from AnalyzeResults import *
from Helpers import *

def get_immediate_subdirectories(dir):
    return [name for name in os.listdir(dir) if os.path.isdir(os.path.join(dir, name))]

sys.path += get_immediate_subdirectories("./")

from inferelator_pipeline import *

exp_data_directory = 'datasets/GeneNetWeaver Data Limit Exp_2'
timeseries_filename = 'insilico_size100_1_dream4_timeseries.tsv'
wildtype_filename = 'insilico_size100_1_wildtype.tsv'
knockout_filename = 'insilico_size100_1_knockouts.tsv'
knockdown_filename = 'insilico_size100_1_knockdowns.tsv'
goldstandard_filename = 'insilico_size100_1_goldstandard.tsv'

# Initialize settings file
settings = {}
settings = ReadConfig(settings)
settings["global"]["working_dir"] = os.getcwd() + '/'
settings["global"]["experiment_name"] = "Data_Limit_Experiment-DFG"

# Read in gold standard network
goldnet = Network()

# Set up output directory
t = datetime.now().strftime("%Y-%m-%d_%H.%M.%S")
settings["global"]["output_dir"] = settings["global"]["output_dir"] + "/" + \
    settings["global"]["experiment_name"] + "-" + t + "/"
os.mkdir(settings["global"]["output_dir"])

# Read in configs for this algorithm
from dfg4grn import *
settings = ReadConfig(settings, "./config/default_values/dfg4grn.cfg")
settings = ReadConfig(settings, settings["dfg4grn"]["config"])

data = {}
knockouts = {}
wildtypes = {}
knockdowns = {}
goldnets = {}
# Loop over the directories we want, reading in the timeseries files
for name in os.listdir(exp_data_directory):
    data[name] = ReadData(exp_data_directory + '/' + name + '/' + timeseries_filename, "timeseries")
    knockouts[name] = ReadData(exp_data_directory + '/' + name + '/' + knockout_filename, "knockout")
    knockdowns[name] = ReadData(exp_data_directory + '/' + name + '/' + knockdown_filename, "knockdown")
    wildtypes[name] = ReadData(exp_data_directory + '/' + name + '/' + wildtype_filename, "wildtype")
    goldnets[name] = exp_data_directory + '/' + name + '/' + goldstandard_filename

jobman = JobManager(settings)

for name in data.keys():
    #dfg = DFG4GRN()
    ts_storage = data[name]
    #settings["dfg4grn"]["eta_z"] = 0.1
    #settings["dfg4grn"]["lambda_w"] = 0.001
    #settings["dfg4grn"]["tau"] = 3
    #settings["global"]["time_series_delta_t"] = (1000 / (len(ts_storage[0].experiments)-1)) / 60.0
    #print settings["global"]["time_series_delta_t"]
    #dfg.setup(ts_storage, TFList(ts_storage[0].gene_list), settings, "DFG-{1}_LambdaW-{0}".format(0.1, name), 20)
    #jobman.queueJob(dfg)

    infjob = InferelatorPipeline()
    infjob.setup(knockouts[name], wildtypes[name], settings, ts_storage, None, "InferelatorPipeline-{0}".format(name))
    jobman.queueJob(infjob)
    infjob = InferelatorPipeline()
    infjob.setup(knockouts[name], wildtypes[name], settings, ts_storage, knockdowns[name], "InferelatorPipeline-{0}-{1}".format(name, "with_kd"))
    jobman.queueJob(infjob)
    infjob = InferelatorPipeline()
    infjob.setup(None, wildtypes[name], settings, ts_storage, None, "InferelatorPipeline-{0}-{1}".format(name, "no_ko"))
    jobman.queueJob(infjob)
    infjob = InferelatorPipeline()
    infjob.setup(knockouts[name], None, settings, ts_storage, None, "InferelatorPipeline-{0}-{1}".format(name, "no_wt"))
    jobman.queueJob(infjob)

jobman.runQueue()
jobman.waitToClear()

goldnet = Network()
goldnet.read_goldstd(goldnets[data.keys()[0]])

tprs, fprs, rocs = GenerateMultiROC(jobman.finished, goldnet, False, settings["global"]["output_dir"] + "/OverallROC.pdf")
ps, rs, precs = GenerateMultiPR(jobman.finished, goldnet, False, settings["global"]["output_dir"] + "/OverallPR.pdf")

SaveResults(jobman.finished, goldnet, settings)

results = []
for job in jobman.finished:
    if job.alg.alg_name == "dfg4grn":
        results.append([job.alg.name, job.alg.best_sign, job.alg.best_sign_all, job.alg.avg_sign, job.alg.avg_sign_all])

out = open(settings["global"]["output_dir"] + "/results.txt",'w')
for r in results:
    out.write(r[0] + "\t" + str(r[1]) + "\t" + str(r[2]) + "\t"+str(r[3])+"\t"+str(r[4])+"\n")

out = open(settings["global"]["output_dir"] + "/results-sorted.txt",'w')
results.sort(key=lambda x: x[0])
for r in results:
    out.write(r[0] + "\t" + str(r[1]) + "\t" + str(r[2]) + "\t"+str(r[3])+"\t"+str(r[4])+"\n")

