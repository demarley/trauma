"""
Created:        16 August 2018
Last Updated:   16 August 2018

Dan Marley
daniel.edison.marley@cernSPAMNOT.ch
Texas A&M University
-----

Base class for plotting deep learning information & performance

Designed for running on desktop at TAMU
with specific set of software installed
--> not guaranteed to work in CMSSW environment!
"""
import os
import sys
import json
import util
import itertools
from datetime import date
import numpy as np

from sklearn.metrics import auc

# load hepPlotter code
try:
    CMSSW_BASE = os.environ['CMSSW_BASE']
    from Analysis.hepPlotter.histogram1D import Histogram1D
    from Analysis.hepPlotter.histogram1D import Histogram2D
    import Analysis.hepPlotter.labels as hpl
    import Analysis.hepPlotter.tools as hpt
except KeyError:
    cwd = os.getcwd()
    hpd = cwd.replace("trauma/miniAna","Analysis/hepPlotter/python/")
    if hpd not in sys.path:
        sys.path.insert(0,hpd)
    from histogram1D import Histogram1D
    from histogram2D import Histogram2D
    import labels as hpl
    import tools as hpt

import ROOT
import plotlabels as plb
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter


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

        self.formatter       = FormatStrFormatter('%g')
        self.betterColors    = hpt.betterColors()['linecolors']
        self.sample_labels   = plb.sample_labels()
        self.variable_labels = plb.variable_labels()

        self.msg_svc      = util.VERBOSE()
        self.output_dir   = ''
        self.image_format = 'pdf'
        self.process_label = ''      # if a single process is used for all training, set this

        self.df = None
        self.targets = []
        self.target_pairs = None

        self.CMSlabelStatus = "Simulation Internal"


    def initialize(self,dataframe,targets={},scaled_inputs=False):
        """
        Set parameters of class to make plots

        @param dataframe      The dataframe that contains the physics information
        @param target_names   The names for different targets to use for making the plots
        @param target_values  The values for the different targets, e.g., [0,1,2,...]
        """
        self.df = dataframe
        self.scaled_inputs = scaled_inputs

        self.listOfFeatures     = [i for i in self.df.keys() if i!='target']
        self.listOfFeaturePairs = list(itertools.combinations(self.listOfFeatures,2))

        for i,(n,v) in enumerate(targets.iteritems()):
            tmp    = Target(n)
            tmp.df = self.df.loc[self.df['target']==v]
            tmp.label = self.sample_labels[n].label
            tmp.color = self.betterColors[i]
            tmp.target_value = v
            self.targets.append(tmp)

        # Create unique combinations of the targets in pairs 
        # (to calculate separation between classes)
        self.target_pairs = list(itertools.combinations(targets.keys(),2))

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
            hist = Histogram1D()

            hist.normed  = True
            hist.stacked = False
            hist.binning = self.variable_labels[feature].binning
            hist.x_label = self.variable_labels[feature].label
            hist.y_label = "A.U." if hist.normed else "Events"
            hist.format  = self.image_format
            hist.saveAs  = self.output_dir+"/hist_"+feature+"_"+self.date
            hist.CMSlabel = 'outer'
            hist.CMSlabelStatus = self.CMSlabelStatus

            hist.ratio.value  = "significance"
            hist.ratio.ylabel = r"S/$\sqrt{\text{B}}$"
            hist.ratio.ymin   = {'ymin':0}

            hist.initialize()

            for target in self.targets:
                kwargs = {"draw_type":"step","edgecolor":target.color,"label":target.label}
                hist.Add(target.df[feature],name=target.name,**kwargs)

            # Add ratio plot
            n_extra_text = 0      # auto-adjust label based on number of extra text args
            for pair in self.target_pairs:
                hist.ratio.Add(numerator=pair[0],denominator=pair[1],draw_type='errorbar')

                name_a = pair[0] 
                name_b = pair[1]

                separation = self.separations[feature]['-'.join(pair)]
                hist.extra_text.Add("Sep({0},{1}) = {2:.2f}".format(name_a,name_b,separation),
                                    coords=[0.03,(0.97-0.08*n_extra_text)])
                n_extra_text+=1

            p = hist.execute()
            hist.savefig()


        # plot the 2D features
        for hi,featurepairs in enumerate(self.listOfFeaturePairs):
            xfeature = featurepairs[0]
            yfeature = featurepairs[1]

            for target in self.targets:
                hist = Histogram2D()

                print xfeature,yfeature
                xbins = self.variable_labels[xfeature].binning
                ybins = self.variable_labels[yfeature].binning
                print xbins,ybins

                hist.colormap = 'default'
                hist.colorbar['title'] = "Events"

                try:
                    hist.binning = [xbins.tolist(),ybins.tolist()]
                except:
                    hist.binning = [xbins,ybins]
                hist.x_label = self.variable_labels[xfeature].label
                hist.y_label = self.variable_labels[yfeature].label
                hist.format  = self.image_format
                hist.saveAs  = self.output_dir+"/hist2d_"+target.name+"_"+xfeature+"-"+yfeature+"_"+self.date
                hist.CMSlabel = 'outer'
                hist.CMSlabelStatus = self.CMSlabelStatus
                hist.logplot['data'] = True
                hist.extra_text.Add(self.sample_labels[target.name].label,coords=[0.03,0.97])
                hist.initialize()

                h,binsx,binsy = np.histogram2d(target.df[xfeature],target.df[yfeature],bins=[xbins,ybins])
                dummyx = []
                dummyy = []
                for x in 0.5*(binsx[:-1]+binsx[1:]):
                    for y in 0.5*(binsy[:-1]+binsy[1:]):
                        dummyx.append(x)
                        dummyy.append(y)

                hist.Add([dummyx,dummyy],weights=h.flatten(),name=target.name)

                p = hist.execute()
                hist.savefig()

        return


    def separation(self):
        """Plot the separations between features of the NN"""
        listOfFeatures = list(self.listOfFeatures) #[self.variable_labels[f].label for f in self.listOfFeatures]
        listOfFeaturePairs = list(self.listOfFeaturePairs)
        featurelabels = [self.variable_labels[f].label for f in self.listOfFeatures]

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
            ax.set_yticks(listOfFeatures)
            ax.set_yticklabels([self.variable_labels[f].label for f in listOfFeatures],fontsize=12)
            ax.set_xticklabels([self.formatter(i) for i in ax.get_xticks()])

            # CMS/COM Energy Label + Signal name
            self.stamp_cms(ax)
            self.stamp_energy(ax)
            ax.text(0.03,0.93,"{0} - {1}".format(target_a,target_b),fontsize=16,
                    ha='left',va='bottom',transform=ax.transAxes)

            plt.savefig("{0}.{1}".format(saveAs,self.image_format))
            plt.close()


            ## Two dimensional separation plot ##
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
            cbar.ax.set_yticklabels([i.get_text().replace("$","") for i in cbar.ax.get_yticklabels()])

            # shift location of ticks to center of the bins
            ax.set_xticks(np.arange(len(self.listOfFeatures))+0.5, minor=False)
            ax.set_yticks(np.arange(len(self.listOfFeatures))+0.5, minor=False)
            ax.set_xticklabels(featurelabels, minor=False, ha='right', rotation=70)
            ax.set_yticklabels(featurelabels, minor=False)

            ## CMS/COM Energy Label + Signal name
            self.stamp_cms(ax)
            self.stamp_energy(ax)

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
            cbar.ax.set_yticklabels([i.get_text().replace("$","") for i in cbar.ax.get_yticklabels()])

            labels = [self.variable_labels[feature].label for feature in corrmat.columns.values]
            # shift location of ticks to center of the bins
            ax.set_xticks(np.arange(len(labels))+0.5, minor=False)
            ax.set_yticks(np.arange(len(labels))+0.5, minor=False)
            ax.set_xticklabels(labels, fontsize=12, minor=False, ha='right', rotation=70)
            ax.set_yticklabels(labels, fontsize=12, minor=False)

            ## CMS/COM Energy Label + Signal name
            self.stamp_cms(ax)
            self.stamp_energy(ax)
            ax.text(0.03,0.93,target.label,fontsize=16,ha='left',va='bottom',transform=ax.transAxes)

            plt.savefig("{0}.{1}".format(saveAs,self.image_format),
                        format=self.image_format,dpi=300,bbox_inches='tight')
            plt.close()

        return



    def prediction(self,train_data={},test_data={},i=0):
        """
        Plot the training and testing predictions. 
        To save on memory, pass this TH1s directly, rather than raw values.
        """
        self.msg_svc.INFO("DL : Plotting DNN prediction. ")
        binning = [bb*0.1 for bb in range(11)]

        # Make a plot for each target value (e.g., QCD prediction to be QCD; Top prediction to be QCD, etc)
        hist = Histogram1D()

        hist.normed  = True  # compare shape differences (likely don't have the same event yield)
        hist.format  = self.image_format
        hist.saveAs  = "{0}/hist_DNN_prediction_{1}".format(self.output_dir,self.date)
        hist.binning = binning
        hist.stacked = False
        hist.x_label = "Prediction"
        hist.y_label = "A.U."
        hist.CMSlabel = 'outer'
        hist.CMSlabelStatus   = self.CMSlabelStatus
        hist.legend['fontsize'] = 18

        hist.ratio.value  = "ratio"
        hist.ratio.ylabel = "Train/Test"

        hist.initialize()

        json_data = {}
        for t,target in enumerate(self.targets):

            target_value = target.target_value  # arrays for multiclassification 

            train_t = train_data[target.name]
            test_t  = test_data[target.name]

            train_kwargs = {"draw_type":"step","edgecolor":target.color,
                            "label":target.label+" Train"}
            test_kwargs  = {"draw_type":"stepfilled","edgecolor":target.color,
                            "color":target.color,"linewidth":0,"alpha":0.5,
                            "label":target.label+" Test"}

            hist.Add(train_t,name=target.name+'_train',**train_kwargs) # Training
            hist.Add(test_t,name=target.name+'_test',**test_kwargs)    # Testing

            hist.ratio.Add(numerator=target.name+'_train',denominator=target.name+'_test')

            ## Save data to JSON file
            if isinstance(train_t,ROOT.TH1):
                d_tr = hpt.hist2list(train_t)
                d_te = hpt.hist2list(test_t)
                json_data[target.name+"_train"] = {"binning":d_tr.bins.tolist(),
                                                   "content":d_tr.content.tolist()}
                json_data[target.name+"_test"]  = {"binning":d_te.bins.tolist(),
                                                   "content":d_te.content.tolist()}
            else:
                json_data[target.name+"_train"] = {"binning":hist.binning,
                                                   "content":train_t}
                json_data[target.name+"_test"]  = {"binning":hist.binning,
                                                   "content":test_t}

        # calculate separation between predictions
        for t,target in enumerate(self.target_pairs):
            data_a = json_data[ target[0]+"_test" ]["content"]
            data_b = json_data[ target[1]+"_test" ]["content"]

            separation = util.getSeparation(data_a,data_b)
            hist.extra_text.Add("Test Sep({0},{1}) = {2:.2f}".format(target[0],target[1],separation),
                                coords=[0.03,(0.97-0.08*t)])

            json_data[ '-'.join(target)+"_test" ] = {"separation":separation}

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

        # Plot ROC curve
        roc_auc = auc(fprs,tprs)
        ax.plot(fprs,tprs,label='AUC = {0:.2f}'.format(roc_auc),lw=2)
        # save ROC curve to CSV file (to plot later)
        csv = [ "{0},{1}".format(fp,tp) for fp,tp in zip(fprs,tprs) ]
        util.to_csv("{0}.csv".format(saveAs),csv)

        ax.set_xlim([0.0, 1.0])
        ax.set_ylim([0.0, 1.5])

        ax.set_xlabel(r'Background',ha='right',va='top',position=(1,0))
        ax.set_ylabel(r'Signal',ha='right',va='bottom',position=(0,1))

        ax.set_xticklabels([self.formatter(i) for i in ax.get_xticks()],fontsize=20)
        ax.set_yticklabels([self.formatter(i) for i in ax.get_yticks()],fontsize=20)

        ## CMS/COM Energy Label
        self.stamp_cms(ax)
        self.stamp_energy(ax)

        if accuracy:
            ax.text(0.97,0.03,r"Accuracy = {0:.2f}$\pm${1:.2f}".format(accuracy['mean'],accuracy['std']),ha='right')

        leg = ax.legend()
        leg.draw_frame(False)

        plt.savefig('{0}.{1}'.format(saveAs,self.image_format))
        plt.close()

        return


    def plot_history(self,history,ax=None,key='loss',index=-1):
        """Draw history of model"""
        try:
            loss     = history.history[key]
            val_loss = history.history.get('val_'+key)
        except:
            loss     = history
            val_loss = None

        x = range(1,len(loss)+1)
        label = key.title()
        if index>=0: label += ' {0}'.format(index)
        ax.plot(x,loss,label=label)
        csv = [ "{0},{1}\n".format(i,j) for i,j in zip(x,loss) ]

        if val_loss is not None:
            label = 'Validation {0}'.format(index) if index>=0 else 'Validation'
            ax.plot(x,val_loss,label=label)
            csv += [ "{0},{1}\n".format(i,j) for i,j in zip(x,val_loss) ]

        return csv


    def history(self,history,kfold=-1):
        """Plot history as a function of epoch for model"""
        self.msg_svc.INFO("DL : Plotting loss as a function of epoch number.")

        for key in ['loss','acc']:
            fig,ax = plt.subplots()

            saveAs   = "{0}/history_{1}_{2}".format(self.output_dir,key,self.date)
            csv      = self.plot_history(history,ax=ax,key=key)
            filename = "{0}.csv".format(saveAs)
            util.to_csv(filename,csv)

            ax.set_xlabel('Epoch',fontsize=22,ha='right',va='top',position=(1,0))
            ax.set_ylabel(key.title(),fontsize=22,ha='right',va='bottom',position=(0,1))

            ax.set_xticklabels([self.formatter(i) for i in ax.get_xticks()],fontsize=20)
            ax.set_yticklabels(['']+[self.formatter(i) for i in ax.get_yticks()[1:-1]]+[''],fontsize=20)

            ## CMS/COM Energy Label
            self.stamp_cms(ax)
            self.stamp_energy(ax)

            leg = ax.legend(loc=1,numpoints=1,fontsize=12,ncol=1,columnspacing=0.3)
            leg.draw_frame(False)

            plt.savefig(self.output_dir+'/{0}_epochs_{1}.{2}'.format(key,self.date,self.image_format),
                        format=self.image_format,bbox_inches='tight',dpi=200)
            plt.close()

        return


    def getSeparations(self):
        """Calculate separations between classes for each feature"""
        self.separations = dict( (k,{}) for k in self.listOfFeatures)
        for featurepairs in self.listOfFeaturePairs:
             self.separations['-'.join(featurepairs)] = {}

        # One dimensional separations
        for target in self.target_pairs:
            target_a = [i for i in self.targets if i.name==target[0]][0]
            target_b = [i for i in self.targets if i.name==target[1]][0]

            saveAs = "{0}/separations1D_{1}-{2}_{3}".format(self.output_dir,target_a.name,target_b.name,self.date)
            fcsv = open("{0}.csv".format(saveAs),"w")

            for feature in self.listOfFeatures:
                if feature=='target': continue

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
            saveAs = "{0}/separations2D_{1}-{2}_{3}".format(self.output_dir,target_a.name,target_b.name,self.date)
            fcsv   = open("{0}.csv".format(saveAs),"w")

            for featurepairs in self.listOfFeaturePairs:
                feature_x = featurepairs[0]
                feature_y = featurepairs[1]

                binning_x = self.variable_labels[feature_x].binning
                binning_y = self.variable_labels[feature_y].binning

                # bin the data to make separation calculation simple
                data_a,_,_ = np.histogram2d(target_a.df[feature_x],target_a.df[feature_y],
                                       bins=[binning_x,binning_y],normed=True)
                data_b,_,_ = np.histogram2d(target_b.df[feature_x],target_b.df[feature_y],
                                       bins=[binning_x,binning_y],normed=True)

                separation = util.getSeparation2D(data_a,data_b)
                self.separations['-'.join(featurepairs)]['-'.join(target)] = separation

                # Save separation info to CSV file
                fcsv.write("{0},{1},{2}".format(feature_x,feature_y,separation))
            fcsv.close() 

        return


    def stamp_energy(self,axis,ha='right',coords=[0.99,1.00],fontsize=16,va='bottom'):
        energy_stamp = hpl.EnergyStamp()
        axis.text(coords[0],coords[1],energy_stamp.text,fontsize=fontsize,ha=ha,va=va,transform=axis.transAxes)

        return

    def stamp_cms(self,axis,ha='left',va='bottom',coords=[0.02,1.00],fontsize=16):
        cms_stamp = hpl.CMSStamp(self.CMSlabelStatus)
        axis.text(coords[0],coords[1],cms_stamp.text,fontsize=fontsize,ha=ha,va=va,transform=axis.transAxes)

        return


## THE END ##
