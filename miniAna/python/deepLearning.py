"""
Created:        11 November  2016
Last Updated:   15 February  2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Base class for performing deep learning 

Designed for running on desktop at TAMU
with specific set of software installed
--> not guaranteed to work in CMSSW environment!

Does not use ROOT directly.
Instead, this is setup to use flat ntuples
that are accessed via uproot.


> UPROOT:     https://github.com/scikit-hep/uproot
> KERAS:      https://keras.io/
> TENSORFLOW: https://www.tensorflow.org/
> PYTORCH:    http://pytorch.org/
> LWTNN:      https://github.com/lwtnn/lwtnn

Expandable: Do 'testing' phase later than training phase
            Diagnostics post-training phase
            Different model (PyTorch)
"""
import json
import util
import datetime

from deepLearningPlotter import DeepLearningPlotter

import ROOT
import uproot
import numpy as np
import pandas as pd

import keras
from keras.models import Sequential,model_from_json,load_model
from keras.layers import Dense, Activation
from keras.callbacks import EarlyStopping
from keras.utils.np_utils import to_categorical
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
from sklearn.multiclass import OneVsRestClassifier


# fix random seed for reproducibility
seed = 2018
np.random.seed(seed)


class DeepLearning(object):
    """Deep Learning base class"""
    def __init__(self):
        self.date = datetime.date.today().strftime('%d%b%Y')

        ## Handling NN objects and data -- set in the class
        self.df  = None          # dataframe containing physics information
        self.fpr = None          # ROC curve: false positive rate
        self.tpr = None          # ROC curve: true positive rate
        self.model = None      # Keras model
        self.accuracy  = {'mean':0,'std':0}   # k-fold accuracies
        self.histories = []           # model history (for ecah k-fold)
        self.train_data = {}          # set later
        self.test_data  = {}          # set later
        self.train_predictions = []   # set later
        self.test_predictions  = []   # set later

        ## NN architecture & parameters -- set by config file
        self.treename   = 'features'    # Name of TTree to access in ROOT file (via uproot)
        self.useLWTNN   = True          # export (& load model from) files for LWTNN
        self.dnn_name   = "dnn"         # name to access in lwtnn ('variables.json')
        self.hep_data   = ""            # Name for loading features (physics data) -- assumes all data in one file
        self.model_name = ""            # Name for saving/loading model
        self.output_dir = 'data/dnn/'   # directory for storing NN data
        self.dnn_method = None          # DNN method applied: classification/regression: ['binary','multi','regression']
        self.runDiagnostics = True      # Make plots pre/post training
        self.verbose_level  = 'INFO'
        self.verbose = False

        self.loss    = 'binary_crossentropy' # preferred for binary classification
        self.init    = 'normal'
        self.nNodes  = []
        self.dropout = None
        self.metrics = ['accuracy']
        self.features   = []
        self.epochs     = 1        
        self.optimizer  = 'adam'
        self.input_dim  = 1                  # len(self.features)
        self.output_dim = 1                  # number of output dimensions (# of categories/# of predictions for regression)
        self.batch_size = 32
        self.activations   = ['elu']         # https://keras.io/activations/
        self.nHiddenLayers = 1
        self.earlystopping = {}              # {'monitor':'loss','min_delta':0.0001,'patience':5,'mode':'auto'}
        self.targets = {"bckg":0,"signal":1}


    def initialize(self):   #,config):
        """Initialize a few parameters after they've been set by user"""
        self.msg_svc       = util.VERBOSE()
        self.msg_svc.level = self.verbose_level
        self.msg_svc.initialize()
        self.verbose = not self.msg_svc.compare(self.verbose_level,"WARNING") # verbose if level is <"WARNING"


        # initialize empty dictionaries, lists
        self.test_data  = {'X':[],'Y':[]}
        self.train_data = {'X':[],'Y':[]}
        self.test_predictions  = []
        self.train_predictions = []

        self.fpr = []  # false positive rate
        self.tpr = []  # true positive rate
        self.histories  = []
        if not self.model_name:  self.model_name = self.hep_data.split('/')[-1].split('.')[0]+'_'+self.date

        ## -- Plotting framework
        print " >> Store output in ",self.output_dir
        self.plotter = DeepLearningPlotter()  # class for plotting relevant NN information
        self.plotter.output_dir   = self.output_dir
        self.plotter.image_format = 'pdf'

        ## -- Adjust model architecture parameters (flexibilty in config file)
        if len(self.nNodes)==1 and self.nHiddenLayers>0:
            # All layers (initial & hidden) have the same number of nodes
            self.msg_svc.DEBUG("DL : Setting all layers ({0}) to have the same number of nodes ({1})".format(self.nHiddenLayers+1,self.nNodes))
            nodes_per_layer = self.nNodes[0]
            self.nNodes = [nodes_per_layer for _ in range(self.nHiddenLayers+1)] # 1st layer + nHiddenLayers

        ## -- Adjust activation function parameter (flexibilty in config file)
        if len(self.activations)==1:
            # Assume the same activation function for all layers (input,hidden,output)
            self.msg_svc.DEBUG("DL : Setting input, hidden, and output layers ({0}) \n".format(self.nHiddenLayers+2)+\
                               "     to have the same activation function {0}".format(self.activations[0]) )
            activation = self.activations[0]
            self.activations = [activation for _ in range(self.nHiddenLayers+2)] # 1st layer + nHiddenLayers + output
        elif len(self.activations)==2 and self.nHiddenLayers>0:
            # Assume the last activation is for the output and the first+hidden layers have the first activation
            self.msg_svc.DEBUG("DL : Setting input and hidden layers ({0}) to the same activation function, {1},\n".format(self.nHiddenLayers+1,self.activations[0])+\
                               "     and the output activation to {0}".format(self.activations[1]) )
            first_hidden_act = self.activations[0]
            output_act       = self.activations[1]
            self.activations = [first_hidden_act for _ in range(self.nHiddenLayers+1)]+[output_act]

        return


    def training(self):
        """Train NN model"""
        self.load_hep_data()
        self.build_model()

        if self.runDiagnostics: self.diagnostics(pre=True)

        self.train_model()      # this also saves the model in case errors occur while plotting evaluations

        # save plots of the performance
        if self.runDiagnostics: self.diagnostics(post=True)

        return


    def inference(self,data=None):
        """
        Run inference of the NN model
        User responsible for diagnostics if not doing training: 
        -> save all predictions (& labels) using 'self.test_predictions'
           then call individual functions:
              plot_features()   -> compare features of the inputs
              plot_prediction() -> compare output prediction (works for classification)
              plot_ROC()        -> signal vs background efficiency (need self.fpr, self.tpr filled)
        """
        self.load_model(self.useLWTNN)

        if data is None:
            try:
                self.load_hep_data()
                data = self.df[self.features]
            except:
                self.msg_svc.ERROR("DL : runInference() cannot proceed because 'data' is None and cannot load HEP data")
                self.msg_svc.ERROR("DL : Please check your implementation.")
                return -999

        prediction = self.predict(data)

        return prediction


    ## Specific functions to perform training/inference tasks
    def build_model(self):
        """Construct the NN model -- only Keras support for now"""
        self.msg_svc.INFO("DL : Build the neural network model")

        ## Declare the model
        self.model = Sequential() # The Keras Sequential model is a linear stack of layers.

        ## Add 1st layer
        self.model.add( Dense( int(self.nNodes[0]), input_dim=self.input_dim, kernel_initializer=self.init, activation=self.activations[0]) )

        ## Add hidden layer(s)
        for h in range(self.nHiddenLayers):
            self.model.add( Dense( int(self.nNodes[h+1]), kernel_initializer=self.init, activation=self.activations[h+1]) )

        ## Add the output layer
        print self.output_dim
        self.model.add( Dense(self.output_dim,kernel_initializer=self.init, activation=self.activations[-1]) )

        ## Build the model
        self.model.compile(loss=self.loss, optimizer=self.optimizer, metrics=self.metrics)

        return



    def train_model(self):
        """Setup for training the model using k-fold cross-validation"""
        self.msg_svc.INFO("DL : Train the model!")

        X_train = self.df[self.features].values  # self.df[self.features].values
        Y_train = self.df['target'].values       # self.df['target'].values

        X_train,X_test,Y_train,Y_test = train_test_split(X_train,Y_train,test_size=0.3) # split into train/test

        print X_train.shape, Y_train.shape

        ## Fit the model to training data & save the history
        callbacks_list = [] if not self.earlystopping else [EarlyStopping(**self.earlystopping)]
        history = self.model.fit(X_train,Y_train,epochs=self.epochs,validation_split=0.25,\
                                 callbacks=callbacks_list,batch_size=self.batch_size,verbose=self.verbose)
        self.histories = history


        # save the model
        self.msg_svc.INFO(" SAVE MODEL")
        self.save_model(self.useLWTNN)


        # evaluate the model
        self.msg_svc.DEBUG("DL :     + Evaluate the model: ")
        train_predictions = self.predict(X_train)

        # Evaluate test sample
        test_predictions  = self.predict(X_test)

        # Make ROC curve from test sample
        fpr,tpr,_ = roc_curve( Y_test, test_predictions )
        self.fpr.append(fpr)
        self.tpr.append(tpr)

        # -- store test/train data from each k-fold as histograms (to compare later)
        h_tests  = dict( (n,ROOT.TH1D("test_"+n,"test_"+n,10,0,1)) for n,v in self.targets.iteritems() )
        h_trains = dict( (n,ROOT.TH1D("train_"+n,"train_"+n,10,0,1)) for n,v in self.targets.iteritems() )

        # fill histogram for each target
        for (n,v) in self.targets.iteritems():
            [h_tests[n].Fill(i)  for i in test_predictions[np.where(Y_test==v)]]
            [h_trains[n].Fill(i) for i in train_predictions[np.where(Y_train==v)]]

        # Plot the predictions to compare test/train
        self.msg_svc.INFO("DL : Plot the train/test predictions")
        self.plotter.prediction(h_trains,h_tests)   # compare DNN prediction for different targets

        self.msg_svc.INFO("DL :   Finished fitting model ")

        return


    def predict(self,data=None):
        """Return the prediction from a test sample"""
        self.msg_svc.DEBUG("DL : Get the DNN prediction")
        if data is None:
            self.msg_svc.ERROR("DL : predict() given NoneType data. Returning -999.")
            self.msg_svc.ERROR("DL : Please check your configuration!")
            return -999.
        return self.model.predict( data )


    def load_hep_data(self,variables2plot=[]):
        """
        Load the physics data (flat ntuple) for NN using uproot
        Convert to DataFrame for easier slicing 

        @param variables2plot    If there are extra variables to plot, 
                                 that aren't features of the NN, include them here
        """
        self.msg_svc.INFO("DL : Load HEP data")
        file = uproot.open(self.hep_data)
        data = file[self.treename]
        df   = data.pandas.df( self.features+['target']+variables2plot )

        self.metadata = file['metadata']   # names of samples, target values, etc.

        # Make the dataset sizes equal (trim away some background)
        fraction = 1
        target_dfs = []
        for k,v in self.targets.iteritems():
            tmp = df[ df.target==v ]
            target_dfs.append(tmp)

        # Find minimum size, shuffle entries
        min_size = min( [k.shape[0] for k in target_dfs] )
        for td,target_df in enumerate(target_dfs):
            # shuffle entries and select first events up to 'min_size'
            if target_df.shape[0]>min_size:
                target_dfs[td] = target_dfs[td].sample(frac=1)[0:min_size]

        self.df = pd.concat( target_dfs )     # re-combine into single dataframe
        self.df = self.df.sample(frac=1)  # shuffle entries

        # set up plotter and save plots of the features and model architecture
        self.plotter.initialize(self.df,self.targets)
        if self.runDiagnostics: self.diagnostics(pre=True)

        return


    def load_model(self,from_lwtnn=True):
        """Load existing model to make plots or predictions"""
        self.model = None

        if from_lwtnn:
            model_json = open(self.model_name+"_model.json",'r').read()
            self.model = model_from_json(model_json)
            self.model.load_weights(self.model_name+"_weights.h5")
            self.model.compile(loss=self.loss, optimizer=self.optimizer, metrics=self.metrics)
        else:
            self.model = load_model('{0}.h5'.format(self.model_name))

        return


    def save_model(self,to_lwtnn=True):
        """Save the model for use later"""
        output = self.output_dir+'/'+self.model_name

        if to_lwtnn:
            ## Save to format for LWTNN
            self.save_features()            ## Save variables to JSON file

            ## model architecture
            model_json = self.model.to_json()
            with open(output+'_model.json', 'w') as outfile:
                outfile.write(model_json)

            ## save the model weights
            self.model.save_weights(output+'_weights.h5')
        else:
            self.model.save('{0}.h5'.format(output))     # creates a HDF5 file of model

        return


    def save_features(self):
        """
        Save the features to a json file to load via lwtnn later
        Hard-coded scale & offset; must change later if necessary
        """
        text = """  {
    "inputs": ["""

        for fe,feature in enumerate(self.features):
            comma = "," if fe!=len(self.features) else ""
            tmp = """
      {"name": "%(feature)s",
       "scale":  1,
       "offset": 0}%(comma)s""" % {'feature':feature,'comma':comma}
            text += tmp
        text += "],"
        text += """
    "class_labels": ["%(name)s"],
    "keras_version": "%(version)s",
    "miscellaneous": {}
  }
""" % {'version':keras.__version__,'name':self.dnn_name}

        varsFileName = self.output_dir+'/variables.json'
        varsFile     = open(varsFileName,'w')
        varsFile.write(text)

        return


    def diagnostics(self,pre=False,post=False):
        """Diagnostic tests of the NN"""

        self.msg_svc.INFO("DL : Diagnostics")

        # Diagnostics before the training
        if pre:
            self.msg_svc.INFO("DL : -- pre-training")
            print ' features '
            self.plotter.features()          # compare features for different targets
            print ' correlations '
            self.plotter.correlation()       # correlations of features
            self.plotter.separation()        # separatiions between classes for features

        # post training/testing
        if post:
            self.msg_svc.INFO("DL : -- post-training :: ROC")
            self.plotter.ROC(self.fpr,self.tpr,self.accuracy)  # ROC curve for signal vs background
            self.msg_svc.INFO("DL : -- post-training :: History")
            self.plotter.loss_history(self.histories) # loss as a function of epoch

        return


## THE END ##
