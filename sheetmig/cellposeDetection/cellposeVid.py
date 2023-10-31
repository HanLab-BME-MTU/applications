from cellpose import models
from cellpose.io import imread
from matplotlib import pyplot as plt
from matplotlib import patches
import os 
import sys
from PyQt5.QtWidgets import QFileDialog, QDialog,QPushButton,QVBoxLayout,QWidget,QApplication,QLabel,QProgressBar,QTextEdit
from PyQt5 import QtGui
import glob
import torch
import pims
import time
import numpy as np
import pandas as pd 
import trackpy as tp
from scipy.spatial import ConvexHull
from tqdm import tqdm
import cv2 as cv
dir = os.path.dirname(os.path.realpath(__file__))+'/'
from matplotlib.colors import hsv_to_rgb

def get_color_from_angle(angle):
    """Get a color corresponding to an angle using HSV colormap"""
    # Convert angle from [-180, 180] to [0, 1]
    hue = (angle + 180) / 360.0
    # Convert hue to RGB color
    return hsv_to_rgb((hue, 1, 1))

def approxEllipseFromArea(img,frameNum):
    ds={"y":0,"x":0,"theta":0,"width":0,"height":0,"ap":0,"thetarads":0,"thetaunshifted":0}
    imgC=np.array(img,dtype=np.uint8)
    if np.max(imgC)!=0:
        cnt,h=cv.findContours(imgC, cv.RETR_EXTERNAL, cv.CHAIN_APPROX_NONE)
        i=np.argmax([len(j) for j in cnt])
        if len(cnt[i]) >= 5:
            el=cv.fitEllipse(cnt[i])
            ds['y']=el[0][0]
            ds['x']=el[0][1]
            ds['theta'] = ((el[2] - 270) % 360) - 180 #rotate theta 90 degrees clockwise such that 0 degrees corresponds to a cell pointing along flow direction (due east)
            #ds['theta']=el[2]
            ds['width']=el[1][0]
            ds['height']=el[1][1]
            apTheta = (el[2] - 90) % 360
            ds['thetaunshifted'] = apTheta
            apRads = np.deg2rad(apTheta)
            ds['thetarads'] = apRads
            ds['ap']= 2 * ((np.cos(apRads) ** 2) - 0.5)    
        else:
            ds['y']=np.nan
            ds['x']=np.nan
            ds['theta']=np.nan
            ds['width']=np.nan
            ds['height']=np.nan
            ds['ap']=np.nan
            print("No ellipse was created for contour. Total contour points: ",len(cnt[i]))
            print("Ellipse creation failed on frame: ",frameNum)
#        if all(el):
#            ds['y']=el[0][0]
#            ds['x']=el[0][1]
#            ds['theta'] = ((el[2] - 270) % 360) - 180 #rotate theta 90 degrees clockwise such that 0 degrees corresponds to a cell pointing along flow direction (due east)
#            #ds['theta']=el[2]
#            ds['width']=el[1][0]
#            ds['height']=el[1][1]
#            apTheta = (el[2] - 90) % 360
#            apRads = np.deg2rad(apTheta)
#            ds['ap']= 2 * (np.cos(apTheta) ** 2 - 0.5)
#        else:
#            ds['y']=np.nan
#            ds['x']=np.nan
#            ds['theta']=np.nan
#            ds['width']=np.nan
#            ds['height']=np.nan
#            ds['ap']=np.nan
#            print("No ellipse was created for contour",cnt[i])
    return ds



class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.resize(800, 600)
        self.model=''
        self.folder=''
        self.mf=0
        self.ff=0
        self.sp=0
        self.consoleText=''
        #file=self.getFolder()
        self.folderButton = QPushButton('Open Project Folder')
        self.folderButton.clicked.connect(self.getFolder)
        self.folderLab=QLabel("No Folder Selected")
        self.modelButton=QPushButton('Select Model File')
        self.modelButton.clicked.connect(self.getModel)
        self.modelLab=QLabel("No Model Selected")
        self.startButton=QPushButton('Start')
        self.startButton.clicked.connect(self.start)
        self.progress=QProgressBar(self)
        self.console=QTextEdit(self.consoleText)
        self.showPlots=QPushButton("Show Plots")
        self.showPlots.clicked.connect(self.plot)
        layout = QVBoxLayout()
        layout.addWidget(self.folderButton)
        layout.addWidget(self.folderLab)
        layout.addWidget(self.modelButton)
        layout.addWidget(self.modelLab)
        layout.addWidget(self.startButton)
        layout.addWidget(self.progress)
        layout.addWidget(self.console)
        layout.addWidget(self.showPlots)

        self.setLayout(layout)

    def log(self,msg):
        self.console.append(msg)
        self.repaint()

    def getFolder(self):
        folder=QFileDialog.getExistingDirectory(self,'Image Directory')
        if folder:
            self.ff=1
            self.folderLab.setText(folder)
            self.folder=folder
        else:
            self.ff=0
            self.folderLab.setText("No Folder Selected")
    def getModel(self):
        [model,_]=QFileDialog.getOpenFileName(self,"Select Model",filter="*")
        if model:
            self.mf=1
            self.modelLab.setText(model)
            self.model=model
        else:
            self.mf=0
            self.modelLab.setText("No Model Selected")
    def start(self):
        if not self.mf or not self.ff:
            print(self.mf)
            print(self.ff)
            print( not self.mf or not self.ff)
            self.log("Model or Filder not selected")
            return
        self.log("Starting...")
        
        #time.sleep(1)
        self.md = models.CellposeModel(pretrained_model=self.model,gpu=True)
        self.log("Loaded Model")
        video=pims.open(self.folder+'/*.tif')
        self.frames=np.array(video)
        self.log("Loaded Frames")
        ch=[0,0]
        self.m=[]
        self.df=pd.DataFrame(columns=["y","x","theta","width","height","ap","thetarads","thetaunshifted","frame"])
        n=len(self.frames)
        self.progress.setMaximum(n)
        self.log("Running model, do not click anything")
        for i,f in enumerate(self.frames):
            masks, flows, styles = self.md.eval(f, diameter=None, channels=ch,progress=True)
            lab=np.unique(masks)
            for l in lab:
                if l!=0:
                    e=approxEllipseFromArea(masks==l,i)
                    e['frame']=i
                    self.df=pd.concat([self.df,pd.DataFrame(e,index=range(1))]).reset_index(drop=True)
            
            self.m.append(masks)
            self.progress.setValue(i)
            self.repaint()
            time.sleep(0.05)
        self.log("Finished")
        self.progress.setValue(n)
        self.sp=1

    def plot(self):
        if not self.sp:
            self.log("Analysis not finished")
            return
        
        fig,ax=plt.subplots(1,2,tight_layout=True,figsize=(10,5),subplot_kw={'projection':'polar'})
        fig2,ax2=plt.subplots(1,2,tight_layout=True,figsize=(10,5))
        theta=[]
        ratio=[]
        # for i,r in self.df.query("frame=={}".format(0)).iterrows():
        #     theta.append(np.round(np.deg2rad(r["theta"]),2))
        #     #ma=max((r["height"],r["width"]))
        #     #mi=min((r["height"],r["width"]))
        #     #ratio.append(ma/mi)
        #     shiftedTheta = r['theta'] + 90 #shift ellipse angle by 90 degrees to realign
        #     e=patches.Ellipse((r['y'],r['x']), r['width'], r['height'],shiftedTheta)
        #     ax2[0].add_artist(e)
        # ax[0].hist(theta,width=2*np.pi/20,bins=20)
        # ax2[0].set_axis_off()
        # ax2[0].imshow(self.frames[0])
        # theta=[]
        # ratio=[]
        # for i,r in self.df.query("frame=={}".format(len(self.frames)-1)).iterrows():
        #     theta.append(np.round(np.deg2rad(r["theta"]),2))
        #     #ma=max((r["height"],r["width"]))
        #     #mi=min((r["height"],r["width"]))
        #     #ratio.append(ma/mi)
        #     shiftedThetaFinal = r['theta'] + 90 #shift ellipse angle by 90 degrees to realign
        #     e=patches.Ellipse((r['y'],r['x']), r['width'], r['height'],shiftedThetaFinal,alpha=0.5)
        #     ax2[1].add_artist(e)
        # ax[1].hist(theta,width=2*np.pi/20,bins=20)
        # ax2[1].set_axis_off()
        # ax2[1].imshow(self.frames[-1])
        for i, r in self.df.query("frame=={}".format(0)).iterrows():
            shiftedTheta = r['theta'] + 90
            color = get_color_from_angle(r['theta'])
            e = patches.Ellipse((r['y'], r['x']), r['width'], r['height'], shiftedTheta, edgecolor=color,
                                facecolor='none', alpha=0.5)
            ax2[0].add_artist(e)
        ax2[0].imshow(self.frames[0], cmap='gray')  # Convert the image to grayscale

        for i, r in self.df.query("frame=={}".format(len(self.frames) - 1)).iterrows():
            shiftedThetaFinal = r['theta'] + 90
            color = get_color_from_angle(r['theta'])
            e = patches.Ellipse((r['y'], r['x']), r['width'], r['height'], shiftedThetaFinal, edgecolor=color,
                                facecolor='none', alpha=0.5)
            ax2[1].add_artist(e)
        ax2[1].imshow(self.frames[-1], cmap='gray')  # Convert the image to grayscale

        fig3,ax3=plt.subplots(1,1,tight_layout=True,figsize=(5,5))
        y=[]
        yerr=[]
        for i in range(len(self.frames)):
            y.append(self.df.query("frame=={}".format(i))['theta'].mean())
            angle_data = self.df.query("frame=={}".format(i))['theta']
            if np.size(angle_data) == 0: #CHECK IF NO CELLS DETECTED TO PREVENT DIVIDE BY 0 IN Standard Error Mean CALCULATION
                sem_mean_angle = np.nan
            else:
                sem_mean_angle = np.std(angle_data,ddof=1) / np.sqrt(np.size(angle_data)) 
            yerr.append(sem_mean_angle)
            #yerr.append(self.df.query("frame=={}".format(i))['theta'].std())
        mean_angle_list = y
        sem_mean_angle_list = yerr
        x=range(0,len(self.frames))
        ax3.errorbar(x,y,yerr=yerr,ecolor='c',elinewidth=0.75,capsize=2)
        ax3.set_xlabel('Frame')
        ax3.set_ylabel('Mean Angle (Deg)')
        ax3.set_ylim(-180,180)
        #fig.savefig(dir+'radial.tif',dpi=800)
        #fig2.savefig(dir+'ellipse.tif',dpi=800)
        #fig3.savefig(dir+'time.tif',dpi=800)

        fig4,ax4=plt.subplots(1,1,tight_layout=True,figsize=(5,5))
        self.tracking=tp.link(self.df,15,memory=10)
        for p in self.tracking["particle"].unique():
            ps=self.tracking.query("particle=={}".format(p))
            ax4.plot(ps["frame"],ps["theta"])
        ax4.set_xlabel('Frame')
        ax4.set_ylabel('Theta')
        ax4.set_ylim(-180,180)

        csv_file_path = 'C:\\Users\\sasevera\\Desktop\\alignmentParamOut.csv'   
        columns_to_write = self.df[['y','x','theta','width','height','thetaunshifted','thetarads','ap','frame']]
        columns_to_write['frame'] = columns_to_write['frame'] + 1
        columns_to_write.to_csv(csv_file_path,index=False)
        fig5,ax5=plt.subplots(1,1,tight_layout=True,figsize=(5,5))
        y=[]
        yerr=[]
        for i in range(len(self.frames)):
            y.append(self.df.query("frame=={}".format(i))['ap'].mean())
            ap_data = self.df.query("frame=={}".format(i))['ap']
            sem_mean_ap = np.std(ap_data,ddof=1) / np.sqrt(np.size(ap_data))
            yerr.append(sem_mean_ap)
            #yerr.append(self.df.query("frame=={}".format(i))['ap'].std())
        frame_list = range(1,int(len(self.frames)+1))
        dfAP = pd.DataFrame({'Frame':frame_list,'Mean Angle':mean_angle_list,'Mean Angle SEM':sem_mean_angle_list,'Mean AP':y,'Mean AP SEM':yerr})
        dfAP.to_csv('C:\\Users\\sasevera\\Desktop\\alignmentParamAveragedOutput.csv',index=False)
        x=range(0,len(self.frames))
        ax5.errorbar(x,y,yerr=yerr,ecolor='c',elinewidth=0.75,capsize=2)
        ax5.set_xlabel('Frame')
        ax5.set_ylabel('Mean Alignment Parameter')
        ax5.set_ylim(-1.1,1.1)
        
        frame_last = len(self.frames)
        fig6,ax6=plt.subplots(1,5,tight_layout=True,figsize=(25,5))
        theta=[]
        #FRAME ONE
        for i,r in self.df.query("frame=={}".format(0)).iterrows():
            theta.append(np.round(np.deg2rad(r["theta"]),2))
            shiftedTheta = r['theta'] + 90 #shift ellipse angle by 90 degrees to realign
            color = get_color_from_angle(r['theta'])
            e = patches.Ellipse((r['y'], r['x']), r['width'], r['height'], shiftedTheta, edgecolor=color,
                                facecolor='none', alpha=0.5)
            ax6[0].add_artist(e)
        ax6[0].set_axis_off()
        ax6[0].imshow(self.frames[0], cmap='gray')
        ax6[0].set_title("Frame #1")
        theta=[]

        #FRAME N*1/5
        frame_n1 = int(np.floor(frame_last * 1 / 4))
        for i,r in self.df.query("frame=={}".format(frame_n1)).iterrows():
            theta.append(np.round(np.deg2rad(r["theta"]),2))
            shiftedThetaFinal = r['theta'] + 90 #shift ellipse angle by 90 degrees to realign
            color = get_color_from_angle(r['theta'])
            e = patches.Ellipse((r['y'], r['x']), r['width'], r['height'], shiftedTheta, edgecolor=color,
                                facecolor='none', alpha=0.5)
            ax6[1].add_artist(e)
        ax6[1].set_axis_off()
        ax6[1].imshow(self.frames[frame_n1], cmap='gray')
        ax6[1].set_title("Frame #{}".format(frame_n1))
        theta=[]
        #FRAME N*2/5
        frame_n2 = int(np.floor(frame_last * 2 / 4))
        for i,r in self.df.query("frame=={}".format(frame_n2)).iterrows():
            theta.append(np.round(np.deg2rad(r["theta"]),2))
            #ma=max((r["height"],r["width"]))
            #mi=min((r["height"],r["width"]))
            shiftedThetaFinal = r['theta'] + 90 #shift ellipse angle by 90 degrees to realign
            color = get_color_from_angle(r['theta'])
            e = patches.Ellipse((r['y'], r['x']), r['width'], r['height'], shiftedTheta, edgecolor=color,
                                facecolor='none', alpha=0.5)
            ax6[2].add_artist(e)
        ax6[2].set_axis_off()
        ax6[2].imshow(self.frames[frame_n2], cmap='gray')
        ax6[2].set_title("Frame #{}".format(frame_n2))
        theta=[]
        #FRAME N*3/4
        frame_n3 = int(np.floor(frame_last * 3 / 4))
        for i,r in self.df.query("frame=={}".format(frame_n3)).iterrows():
            theta.append(np.round(np.deg2rad(r["theta"]),2))
            #ma=max((r["height"],r["width"]))
            #mi=min((r["height"],r["width"]))
            shiftedThetaFinal = r['theta'] + 90 #shift ellipse angle by 90 degrees to realign
            color = get_color_from_angle(r['theta'])
            e = patches.Ellipse((r['y'], r['x']), r['width'], r['height'], shiftedTheta, edgecolor=color,
                                facecolor='none', alpha=0.5)
            ax6[3].add_artist(e)
        ax6[3].set_axis_off()
        ax6[3].imshow(self.frames[frame_n3], cmap='gray')
        ax6[3].set_title("Frame #{}".format(frame_n3))
        theta=[]
        #FRAME N-1
        for i,r in self.df.query("frame=={}".format(frame_last-1)).iterrows():
            theta.append(np.round(np.deg2rad(r["theta"]),2))
            #ma=max((r["height"],r["width"]))
            #mi=min((r["height"],r["width"]))
            shiftedThetaFinal = r['theta'] + 90 #shift ellipse angle by 90 degrees to realign
            color = get_color_from_angle(r['theta'])
            e = patches.Ellipse((r['y'], r['x']), r['width'], r['height'], shiftedTheta, edgecolor=color,
                                facecolor='none', alpha=0.5)
            ax6[4].add_artist(e)
        ax6[4].set_axis_off()
        ax6[4].imshow(self.frames[-1], cmap='gray')
        ax6[4].set_title("Frame #{}".format(frame_last))

        # Create an axis for the colorbar next to ax6
        angle = (180 / np.pi) * np.angle(h)
        cbar_ax = fig6.add_axes([0.92, 0.15, 0.02, 0.25])  # [left, bottom, width, height]

        # Create the colorbar
        cbar = fig6.colorbar(im, cax=cbar_ax, orientation='vertical', ticks=[angle.min(), angle.max()])
        cbar.set_label('Angle (degrees)')

        # Display the min and max values on the colorbar
        cbar.ax.set_yticklabels(
            [f'{angle.min():.2f}', f'{angle.max():.2f}'])  # Display with two decimal places

        plt.show()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    files=[]
    main = MainWindow()
    main.show()
    sys.exit(app.exec())