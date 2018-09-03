"""
Created:        11 November  2016
Last Updated:   16 February  2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Base class for plotting deep learning information & performance

Designed for running on desktop at TAMU
with specific set of software installed
--> not guaranteed to work in CMSSW environment!
"""
import json
import util
import itertools
from datetime import date
import numpy as np

from Analysis.hepPlotter.histogram1D import Histogram1D
from Analysis.hepPlotter.histogram1D import Histogram2D

import Analysis.hepPlotter.labels as hpl
import Analysis.hepPlotter.tools as hpt



class Target(object):
    """Class to contain information for targets used in training/inference"""
    def __init__(self,name=""):
        self.name  = name     # Name of this target, e.g., 'signal'
        self.df    = None     # dataframe of this target's features
        self.color = 'k'
        self.label = ''
        self.target_value = -999
        self.binning = 1



class DeepLearningPlotter(object):
    """Plotting utilities for deep learning"""
    def __init__(self):
        """Give default values to member variables"""
        self.date = date.today().strftime('%d%b%Y')

        self.betterColors    = hpt.betterColors()['linecolors']
        self.sample_labels   = hpl.sample_labels()
        self.variable_labels = hpl.variable_labels()

        self.msg_svc      = util.VERBOSE()
        self.output_dir   = ''
        self.image_format = 'pdf'
        self.process_label = ''      # if a single process is used for all training, set this

        self.df = None
        self.targets = []
        self.target_pairs = None

        self.CMSlabelStatus = "Simulation Internal"


    def initialize(self,dataframe,target_names=[],target_values=[]):
        """
        Set parameters of class to make plots

        @param dataframe      The dataframe that contains the physics information
        @param target_names   The names for different targets to use for making the plots
        @param target_values  The values for the different targets, e.g., [0,1,2,...]
        """
        self.df = dataframe

        self.listOfFeatures     = self.df.keys()
        self.listOfFeaturePairs = list(itertools.combinations(self.listOfFeatures,2))

        assert len(target_names)==len(target_values), "Lengths of target names and target values are not the same"

        for i,(n,v) in enumerate(zip(target_names,target_values)):
            tmp    = Target(n)
            tmp.df = self.df.loc[self.df['target']==v]
            tmp.label = self.sample_labels[n].label
            tmp.color = self.betterColors[i]
            tmp.target_value = v
            self.targets.append(tmp)

        # Create unique combinations of the targets in pairs 
        # (to calculate separation between classes)
        target_names = [i.name for i in self.targets]
        self.target_pairs = list(itertools.combinations(target_names,2))

        self.getSeparations()

        return


    def features(self):
        """
        Plot the features
          For classification, compare different targets
          For regression, just plot the features        <- should do data/mc plots instead!
        """
        self.msg_svc.INFO("DL : Plotting features.")

        # plot the features and calculate significance
        for hi,feature in enumerate(self.listOfFeatures):

            hist = HepPlotterHist1D()

            hist.normed  = True
            hist.stacked = False
            hist.binning = self.variable_labels[feature].binning
            hist.x_label = self.variable_labels[feature].label
            hist.y_label = "A.U." if hist.normed else "Events"
            hist.format  = self.image_format
            hist.saveAs  = self.output_dir+"/hist_"+feature+"_"+self.date
            hist.ratio_plot     = True
            hist.ratio_type     = 'significance'
            hist.y_label_ratio  = r"S/$\sqrt{\text{B}}$"
            hist.CMSlabel       = 'top left'
            hist.CMSlabelStatus = self.CMSlabelStatus

            hist.initialize()

            # Add some extra text to the plot
            n_extra_text = 0      # auto-adjust label based on number of extra text args
            if self.processlabel: # physics process that produces these features
                hist.extra_text.Add(self.processlabel,coords=[0.03,0.80])
                n_extra_text+=1


            # Plot the distribution for each target with ratios between classes
            for t,target in enumerate(self.targets):
                ratios = []
                name = target.name
                for pair in self.target_pairs:
                    # only check with first entry to prevent double-plotting ratios
                    if name == pair[0]: ratios.append( (pair[1],True))

                kwargs = {"draw_type":"step","edgecolor":target.color,"label":target.label}

                hist.Add(target.df[feature],name=name,ratios=ratios,**kwargs)


            # Add text to display the separation between each class for the feature
            for target in self.target_pairs:
                name_a = target[0] 
                name_b = target[1]

                separation = self.separations[feature]['-'.join(target)]

                hist.extra_text.Add("Sep({0},{1}) = {2:.2f}".format(name_a,name_b,separation),
                                    coords=[0.03,(0.80-0.08*n_extra_text)])
                n_extra_text+=1

            p = hist.execute()
            hist.savefig()

        return


    def separation(self):
        """Plot the separations between features of the NN"""
        listOfFeatures = list(self.listOfFeatures)
        listOfFeaturePairs = list(self.listOfFeaturePairs)

        nfeatures = len(listOfFeatures)

        for target in self.target_pairs:
            target_a = target[0]
            target_b = target[1]

            ## One dimensional separation plot (horizontal bar chart)
            saveAs = "{0}/separations1D_{1}-{2}_{3}".format(self.output_dir,target_a,target_b,self.date)

            separations = [self.separations[f]['-'.join(target)] for f in listOfFeatures]

            # sort data by separation value
            data = list( zip(listOfFeatures,separations) )
            data.sort(key=lambda x: x[1])
            listOfFeatures[:],separations[:] = zip(*data)

            # make the bar plot
            fig,ax = plt.subplots()
            ax.barh(listOfFeatures, separations, align='center')

            # CMS/COM Energy Label + Signal name
            self.stamp_cms(ax)
            self.stamp_energy(ax)
            ax.text(0.03,0.93,"{0} - {1}".format(target_a,target_b),fontsize=16,
                    ha='left',va='bottom',transform=ax.transAxes)

            plt.savefig("{0}.{1}".format(saveAs,self.image_format))
            plt.close()


            ## Two dimensional separation plot
            saveAs = "{0}/separations2D_{1}-{2}_{3}".format(self.output_dir,target_a,target_b,self.date)

            # from the separations values for each unique (feature_x,feature_y) combination
            # build a matrix that can be drawn using hist2d()
            separations = []
            x_coord = []
            y_coord = []

            for f in listOfFeaturePairs:
                separations.append(self.separations['-'.join(f)]['-'.join(target)])
                x_coord.append(self.listOfFeatures.index(f[0]))
                y_coord.append(self.listOfFeatures.index(f[1]))

            # Now repeat the entries with flipped indices to get the full matrix
            x = list(x_coord)+list(y_coord)
            y = list(y_coord)+list(x_coord)
            separations += separations

            # Plot separation matrix
            # -- Use matplotlib directly
            fig,ax = plt.subplots()

            plt.hist2d(x,y,bins=[range(nfeatures+1),range(nfeatures+1)],
                       weights=separations,vmin=0.0)
            cbar = plt.colorbar()
            cbar.ax.set_ylabel("Separation")

            # shift location of ticks to center of the bins
            ax.set_xticks(np.arange(len(self.listOfFeatures))+0.5, minor=False)
            ax.set_yticks(np.arange(len(self.listOfFeatures))+0.5, minor=False)
            ax.set_xticklabels(self.listOfFeatures, minor=False, ha='right', rotation=70)
            ax.set_yticklabels(self.listOfFeatures, minor=False)

            ## CMS/COM Energy Label + Signal name
            self.stamp_cms(ax)
            self.stamp_energy(ax)
            ax.text(0.03,0.93,target.label,fontsize=16,ha='left',va='bottom',transform=ax.transAxes)

            plt.savefig("{0}.{1}".format(saveAs,self.image_format))
            plt.close()

        return


    def correlation(self):
        """Plot correlations between features of the NN"""
        opts = {'cmap':plt.get_cmap("bwr"),'vmin':-1,'vmax':1}

        for c,target in enumerate(self.targets):
            saveAs = "{0}/correlations_{1}_{2}".format(self.output_dir,target.name,self.date)

            t_ = target.df.drop('target',axis=1)
            corrmat = t_.corr()

            # Save correlation matrix to CSV file
            corrmat.to_csv("{0}.csv".format(saveAs))

            # Plot correlation matrix
            # -- Use matplotlib directly
            fig,ax = plt.subplots()

            heatmap1 = ax.pcolor(corrmat, **opts)
            cbar     = plt.colorbar(heatmap1, ax=ax)

            labels = [self.variable_labels[feature].label for feature in corrmat.columns.values]
            # shift location of ticks to center of the bins
            ax.set_xticks(np.arange(len(labels))+0.5, minor=False)
            ax.set_yticks(np.arange(len(labels))+0.5, minor=False)
            ax.set_xticklabels(labels, fontProperties, fontsize=12, minor=False, ha='right', rotation=70)
            ax.set_yticklabels(labels, fontProperties, fontsize=12, minor=False)

            ## CMS/COM Energy Label + Signal name
            self.stamp_cms(ax)
            self.stamp_energy(ax)
            ax.text(0.03,0.93,target.label,fontsize=16,ha='left',va='bottom',transform=ax.transAxes)

            plt.savefig("{0}.{1}".format(saveAs,self.image_format),
                        format=self.image_format,dpi=300,bbox_inches='tight')
            plt.close()

        return


    def prediction(self,train_data={},test_data={}):
        """Plot the training and testing predictions"""
        self.msg_svc.INFO("DL : Plotting DNN prediction. ")

        values = {0:np.array([1,0,0]), 1:np.array([0,1,0]), 2:np.array([0,0,1])}

        # Plot all k-fold cross-validation results
        for i,(train,trainY,test,testY) in enumerate(zip(train_data['X'],train_data['Y'],test_data['X'],test_data['Y'])):

            # Make a plot for each target value (e.g., QCD prediction to be QCD; Top prediction to be QCD, etc)
            for t_ind,tar in enumerate(self.targets):

                hist = HepPlotter("histogram",1)

                hist.ratio_plot = "ratio"
                hist.normed  = True  # compare shape differences (likely don't have the same event yield)
                hist.format  = self.image_format
                hist.saveAs  = "{0}/hist_DNN_prediction_kfold{1}_target{2}_{3}".format(self.output_dir,i,t_ind,self.date)
                hist.binning = [bb*0.1 for bb in range(11)]
                hist.stacked = False
                hist.x_label = "Prediction"
                hist.y_label = "A.U."
                hist.y_label_ratio = "Test/Train"
                hist.CMSlabel = 'top left'
                hist.CMSlabelStatus   = self.CMSlabelStatus

                hist.extra_text.Add("{0} Prediction".format(tar.name),coords=[0.03,0.80],fontsize=14)
                if self.processlabel: hist.extra_text.Add(self.processlabel,coords=[0.03,0.72],fontsize=14)

                hist.initialize()

                json_data = {}
                for t,target in enumerate(self.targets):

                    target_value = values[target.target_value]  # arrays for multiclassification 

                    train_t = train[np.where((trainY==target_value).all(axis=1))][:,t_ind] # get the t_ind prediction for this target 
                    test_t  = test[np.where((testY==target_value).all(axis=1))][:,t_ind]

                    train_kwargs = {"draw_type":"step","edgecolor":target.color,
                                    "label":target.label+" Train"}
                    test_kwargs  = {"draw_type":"stepfilled","edgecolor":target.color,
                                    "color":target.color,"linewidth":0,"alpha":0.5,
                                    "label":target.label+" Test"}

                    ## Training
                    hist.Add(train_t,name=target.name+'_train',ratios=[],**train_kwargs)

                    ## Testing
                    ratios = [(target.name+'_train',True)]
                    hist.Add(test_t,name=target.name+'_test',ratios=ratios,**test_kwargs)

                    ## Save data to JSON file
                    json_data[target.name+"_train"] = {}
                    json_data[target.name+"_test"]  = {}
                    d_tr,b_tr = np.histogram(train_t,bins=hist.binning,normed=hist.normed)
                    d_te,b_te = np.histogram(test_t,bins=hist.binning,normed=hist.normed)

                    json_data[target.name+"_train"]["binning"] = b_tr.tolist()
                    json_data[target.name+"_train"]["content"] = d_tr.tolist()
                    json_data[target.name+"_test"]["binning"] = b_te.tolist()
                    json_data[target.name+"_test"]["content"] = d_te.tolist()

                # calculate separation between predictions
                for t,target in enumerate(self.target_pairs):
                    data_a = json_data[ target[0]+"_test" ]["content"]
                    data_b = json_data[ target[1]+"_test" ]["content"]

                    separation = util.getSeparation(data_a,data_b)

                    hist.extra_text.Add("Test Sep({0},{1}) = {2:.2f}".format(target[0],target[1],separation),
                                        coords=[0.03,(0.80-0.08*t)])

                    json_data[ '-'.join(target)+"_test" ]["separation"] = separation

                p = hist.execute()
                hist.savefig()

                # save results to JSON file (just histogram values & bins) to re-make plots
                with open("{0}.json".format(hist.saveAs), 'w') as outfile:
                    json.dump(json_data, outfile)

        return



    def ROC(self,fprs=[],tprs=[],accuracy={}):
        """Plot the ROC curve & save to text file"""
        self.msg_svc.INFO("DL : Plotting ROC curve.")

        saveAs = "{0}/roc_curve_{1}".format(self.output_dir,self.date)

        ## Use matplotlib directly
        fig,ax = plt.subplots()

        # Draw all of the ROC curves from the K-fold cross-validation
        ax.plot([0,1],[0,1],ls='--',label='No Discrimination',lw=2,c='gray')
        ax.axhline(y=1,lw=1,c='k',ls='-')

        for ft,(fpr,tpr) in enumerate(zip(fprs,tprs)):
            # Plot ROC curve
            roc_auc = auc(fpr,tpr)
            ax.plot(fpr,tpr,label='K-fold {0} (AUC = {1:.2f})'.format(ft,roc_auc),lw=2)

            # save ROC curve to CSV file (to plot later)
            outfile_name = "{0}_{1}.csv".format(saveAs,ft)
            csv = [ "{0},{1}".format(fp,tp) for fp,tp in zip(fpr,tpr) ]
            util.to_csv(outfile_name,csv)

        ax.set_xlim([0.0, 1.0])
        ax.set_ylim([0.0, 1.5])

        ax.set_xlabel(r'$\epsilon$(anti-top)',ha='right',va='top',position=(1,0))
        ax.set_ylabel(r'$\epsilon$(top)',ha='right',va='bottom',position=(0,1))

        ## CMS/COM Energy Label
        cms_stamp = hpl.CMSStamp(self.CMSlabelStatus)
        cms_stamp.coords = [0.03,0.97]
        cms_stamp.fontsize = 16
        ax.text(cms_stamp.coords[0],cms_stamp.coords[1],cms_stamp.text,fontsize=cms_stamp.fontsize,
                ha=cms_stamp.ha,va=cms_stamp.va,transform=ax.transAxes)

        energy_stamp    = hpl.EnergyStamp()
        energy_stamp.coords = [0.03,0.90]
        energy_stamp.fontsize = 16
        ax.text(energy_stamp.coords[0],energy_stamp.coords[1],energy_stamp.text, 
                fontsize=energy_stamp.fontsize,ha=energy_stamp.ha, va=energy_stamp.va, transform=ax.transAxes)

        if self.processlabel: 
            ax.text(0.03,0.82,self.processlabel)
        if accuracy:
            ax.text(0.03,0.75,r"Accuracy = {0:.2f}$\pm${1:.2f}".format(accuracy['mean'],accuracy['std']))

        leg = ax.legend()
        leg.draw_frame(False)

        plt.savefig('{0}.{1}'.format(saveAs,self.image_format))
        plt.close()

        return


    def plot_loss_history(self,history,ax=None,index=-1):
        """Draw history of model"""
        loss  = history.history['loss']
        x     = range(1,len(loss)+1)
        label = 'Loss {0}'.format(index) if index>=0 else 'Loss'
        ax.plot(x,loss,label=label)

        csv = [ "{0},{1}".format(i,j) for i,j in zip(x,loss) ]

        return csv


    def loss_history(self,history,kfold=0,val_loss=0.0):
        """Plot loss as a function of epoch for model"""
        self.msg_svc.INFO("DL : Plotting loss as a function of epoch number.")

        saveAs = "{0}/loss_epochs_{1}".format(self.output_dir,self.date)
        all_histories = type(history)==list

        # draw the loss curve
        fig,ax = plt.subplots()

        # also save the data to a CSV file
        if all_histories:
            for i,h in enumerate(history):
                csv = self.plot_loss_history(h,ax=ax,index=i)
                filename = "{0}_{1}.csv".format(saveAs,i)
                util.to_csv(filename,csv)
        else:
            csv = self.plot_loss_history(history,ax=ax)
            filename = "{0}.csv".format(saveAs)
            util.to_csv(filename,csv)

        ax.set_xlabel('Epoch',fontsize=22,ha='right',va='top',position=(1,0))
        ax.set_xticklabels(["{0:.1f}".format(i) for i in ax.get_xticks()],fontsize=22)
        ax.set_ylabel('Loss',fontsize=22,ha='right',va='bottom',position=(0,1))
        ax.set_yticklabels(['']+["{0:.1f}".format(i) for i in ax.get_yticks()[1:-1]]+[''],fontsize=22)


        ## CMS/COM Energy Label
        cms_stamp = hpl.CMSStamp(self.CMSlabelStatus)
        cms_stamp.coords = [0.03,0.97]
        cms_stamp.fontsize = 18
        ax.text(cms_stamp.coords[0],cms_stamp.coords[1],cms_stamp.text,fontsize=cms_stamp.fontsize,
                ha=cms_stamp.ha,va=cms_stamp.va,transform=ax.transAxes)

        energy_stamp    = hpl.EnergyStamp()
        energy_stamp.coords = [0.03,0.90]
        energy_stamp.fontsize = 18
        ax.text(energy_stamp.coords[0],energy_stamp.coords[1],energy_stamp.text, 
                fontsize=energy_stamp.fontsize,ha=energy_stamp.ha, va=energy_stamp.va, transform=ax.transAxes)

        text_args = {'ha':'left','va':'top','fontsize':18,'transform':ax.transAxes}
        text = "Validation Loss = {0}; {1} K-folds".format(val_loss,len(history)) if all_histories else "Validation Loss = {0}".format(val_loss)
        ax.text(0.03,0.76,text,**text_args)

        leg = ax.legend(loc=1,numpoints=1,fontsize=12,ncol=1,columnspacing=0.3)
        leg.draw_frame(False)

        plt.savefig(self.output_dir+'/loss_epochs_{0}_{1}.{2}'.format(kfold,self.date,self.image_format),
                    format=self.image_format,bbox_inches='tight',dpi=200)
        plt.close()


        return


    def getSeparations(self):
        """Calculate separations between classes for each feature"""
        self.separations = {}

        # One dimensional separations
        for target in self.target_pairs:
            target_a = [i for i in self.targets if i.name==target[0]][0]
            target_b = [i for i in self.targets if i.name==target[1]][0]

            saveAs = "{0}/separations2D_{1}-{2}_{3}".format(self.output_dir,target_a,target_b,self.date)
            fcsv = open("{0}.csv".format(saveAs),"w")

            for feature in self.listOfFeatures:

                # bin the data to make separation calculation simple
                data_a,_ = np.histogram(target_a.df[feature],normed=True,
                                        bins=self.variable_labels[feature].binning)
                data_b,_ = np.histogram(target_b.df[feature],normed=True,
                                        bins=self.variable_labels[feature].binning)

                separation = util.getSeparation(data_a,data_b)
                self.separations[feature]['-'.join(target)] = separation

                # Save separation info to CSV file
                fcsv.write("{0},{1}".format(feature,separation))
            fcsv.close() 


            # Two dimensional separations
            saveAs = "{0}/separations2D_{1}-{2}_{3}".format(self.output_dir,target_a,target_b,self.date)
            fcsv   = open("{0}.csv".format(saveAs),"w")

            for featurepairs in self.listOfFeaturePairs:
                feature_x = featurepairs[0]
                feature_y = featurepairs[1]

                binning_x = self.variable_labels[feature_x].binning
                binning_y = self.variable_labels[feature_y].binning

                # bin the data to make separation calculation simple
                data_a,_,_ = np.hist2d(target_a.df[feature_x],target_a.df[feature_y],
                                       bins=[binning_x,binning_y],normed=True)
                data_b,_,_ = np.hist2d(target_b.df[feature_x],target_b.df[feature_y],
                                       bins=[binning_x,binning_y],normed=True)

                separation = util.getSeparation2D(data_a,data_b)
                self.separations['-'.join(featurepairs)]['-'.join(target)] = separation

                # Save separation info to CSV file
                fcsv.write("{0},{1},{2}".format(feature_x,feature_y,separation))
            fcsv.close() 

        return



    def model(self,model,name):
        """Plot the model architecture to view later"""
        keras_plot(model,to_file='{0}/{1}_model.eps'.format(self.output_dir,name),show_shapes=True)
        return


    def stamp_energy(self,axis):
        energy_stamp    = hpl.EnergyStamp()
        energy_stamp.ha = 'right'
        energy_stamp.coords = [0.99,1.00]
        energy_stamp.fontsize = 16
        energy_stamp.va = 'bottom'
        ax.text(energy_stamp.coords[0],energy_stamp.coords[1],energy_stamp.text, 
                fontsize=energy_stamp.fontsize,ha=energy_stamp.ha,va=energy_stamp.va,
                transform=axis.transAxes)

        return

    def stamp_cms(self,axis):
        cms_stamp = hpl.CMSStamp(self.CMSlabelStatus)
        cms_stamp.coords = [0.02,1.00]
        cms_stamp.fontsize = 16
        cms_stamp.va = 'bottom'
        ax.text(0.02,1.00,cms_stamp.text,fontsize=cms_stamp.fontsize,
                ha=cms_stamp.ha,va=cms_stamp.va,transform=axis.transAxes)

        return


## THE END ##
