from cellpose import models
from cellpose.io import imread
from matplotlib import pyplot as plt
from matplotlib import patches
import matplotlib as mpl
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
import pickle
# import dill
dir = os.path.dirname(os.path.realpath(__file__))+'/'
from matplotlib.colors import hsv_to_rgb

def get_color_from_angle(angle):
    """Get a color corresponding to an angle using HSV colormap"""

    # Check if angle is NaN or not a number
    if angle is None or (isinstance(angle, (float, int)) and np.isnan(angle)):
        return (0, 0, 0)  # Return black color for NaN or None

    # Convert angle from [0, 180] to [0, 1]
    hue = (angle % 180) / 180.0

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
            # ds['theta'] = ((el[2] - 270) % 360) - 180 #rotate theta 90 degrees clockwise such that 0 degrees corresponds to a cell pointing along flow direction (due east)
            ds['theta'] = ((el[2] + 90) % 180)  #rotate theta 90 degrees clockwise such that 0 degrees corresponds to a cell pointing along flow direction (due east)
            #ds['theta']=el[2]
            ds['width']=el[1][0]
            ds['height']=el[1][1]
            apTheta = (el[2] - 90) % 180
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
        fig6,ax6=plt.subplots(1,6,tight_layout=True,figsize=(24,4))
        theta=[]
        # Create a list to store the apTheta values for each frame
        all_apTheta = []
        #FRAME ONE
        for i,r in self.df.query("frame=={}".format(0)).iterrows():
            theta.append(np.round(np.deg2rad(r["theta"]),2))
            shiftedTheta = r['theta'] + 90 #shift ellipse angle by 90 degrees to realign
            color = get_color_from_angle(r['theta'])
            e = patches.Ellipse((r['y'], r['x']), r['width'], r['height'], shiftedThetaFinal, edgecolor=color,
                                facecolor='none', alpha=0.5)
            ax6[0].add_artist(e)
            # Append the apTheta values to the list
            all_apTheta.append(r['theta'])
        ax6[0].set_axis_off()
        # Compute the mean and standard deviation for the self.frames[0] matrix
        mu = np.mean(self.frames[0])
        sigma = np.std(self.frames[0])
        # Set a constant factor for contrast adjustment
        c = 1.5
        vmin = mu - c * 0.8 * sigma
        vmax = mu + c * sigma
        # Display the self.frames[0] matrix with the determined contrast values
        ax6[0].imshow(self.frames[0], cmap='gray', aspect='auto', interpolation='none', vmin=vmin, vmax=vmax)

        ax6[0].set_title("Frame #1")
        theta=[]

        #FRAME N*1/5
        frame_n1 = int(np.floor(frame_last * 1 / 4))
        for i,r in self.df.query("frame=={}".format(frame_n1)).iterrows():
            theta.append(np.round(np.deg2rad(r["theta"]),2))
            shiftedThetaFinal = r['theta'] + 90 #shift ellipse angle by 90 degrees to realign
            color = get_color_from_angle(r['theta'])
            e = patches.Ellipse((r['y'], r['x']), r['width'], r['height'], shiftedThetaFinal, edgecolor=color,
                                facecolor='none', alpha=0.5)
            ax6[1].add_artist(e)
            # Append the apTheta values to the list
            all_apTheta.append(r['theta'])
        ax6[1].set_axis_off()
        # Compute the mean and standard deviation for the self.frames[0] matrix
        mu = np.mean(self.frames[frame_n1])
        sigma = np.std(self.frames[frame_n1])
        # Set a constant factor for contrast adjustment
        vmin = mu - c * 0.8 * sigma
        vmax = mu + c * sigma
        # Display the self.frames[0] matrix with the determined contrast values
        ax6[1].imshow(self.frames[frame_n1], cmap='gray', aspect='auto', interpolation='none', vmin=vmin, vmax=vmax)
        # ax6[1].imshow(self.frames[frame_n1], cmap='gray')
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
            e = patches.Ellipse((r['y'], r['x']), r['width'], r['height'], shiftedThetaFinal, edgecolor=color,
                                facecolor='none', alpha=0.5)
            ax6[2].add_artist(e)
            # Append the apTheta values to the list
            all_apTheta.append(r['theta'])
        ax6[2].set_axis_off()
        # Compute the mean and standard deviation for the self.frames[0] matrix
        mu = np.mean(self.frames[frame_n2])
        sigma = np.std(self.frames[frame_n2])
        # Set a constant factor for contrast adjustment
        vmin = mu - c * 0.8 * sigma
        vmax = mu + c * sigma
        # Display the self.frames[0] matrix with the determined contrast values
        ax6[2].imshow(self.frames[frame_n2], cmap='gray', aspect='auto', interpolation='none', vmin=vmin, vmax=vmax)

        # ax6[2].imshow(self.frames[frame_n2], cmap='gray')
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
            e = patches.Ellipse((r['y'], r['x']), r['width'], r['height'], shiftedThetaFinal, edgecolor=color,
                                facecolor='none', alpha=0.5)
            ax6[3].add_artist(e)
            # Append the apTheta values to the list
            all_apTheta.append(r['theta'])
        ax6[3].set_axis_off()
        # Compute the mean and standard deviation for the self.frames[0] matrix
        mu = np.mean(self.frames[frame_n3])
        sigma = np.std(self.frames[frame_n3])
        # Set a constant factor for contrast adjustment
        vmin = mu - c * 0.8 * sigma
        vmax = mu + c * sigma
        # Display the self.frames[0] matrix with the determined contrast values
        ax6[3].imshow(self.frames[frame_n3], cmap='gray', aspect='auto', interpolation='none', vmin=vmin, vmax=vmax)

        # ax6[3].imshow(self.frames[frame_n3], cmap='gray')
        ax6[3].set_title("Frame #{}".format(frame_n3))
        theta=[]

        # try:
        #     with open('self.pkl', 'wb') as f:
        #         pickle.dump(self, f)
        # except Exception as e:
        #     print(f"Error during pickling: {e}")
        #
        # with open('self_dill.pkl', 'wb') as f2:
        #     dill.dump(self, f2)

        #FRAME N-1
        for i,r in self.df.query("frame=={}".format(frame_last-1)).iterrows():
            theta.append(np.round(np.deg2rad(r["theta"]),2))
            #ma=max((r["height"],r["width"]))
            #mi=min((r["height"],r["width"]))
            shiftedThetaFinal = r['theta'] + 90 #shift ellipse angle by 90 degrees to realign
            color = get_color_from_angle(r['theta'])
            e = patches.Ellipse((r['y'], r['x']), r['width'], r['height'], shiftedThetaFinal, edgecolor=color,
                                facecolor='none', alpha=0.5)
            ax6[4].add_artist(e)
            # Append the apTheta values to the list
            all_apTheta.append(r['theta'])
        ax6[4].set_axis_off()
        # Compute the mean and standard deviation for the self.frames[0] matrix
        mu = np.mean(self.frames[-1])
        sigma = np.std(self.frames[-1])
        # Set a constant factor for contrast adjustment
        vmin = mu - c * 0.8 * sigma
        vmax = mu + c * sigma
        # Display the self.frames[0] matrix with the determined contrast values
        im = ax6[4].imshow(self.frames[-1], cmap='gray', aspect='auto', interpolation='none', vmin=vmin, vmax=vmax)
        ax6[4].set_title("Frame #{}".format(frame_last))

        # Create a colormap
        cmap = plt.get_cmap("hsv")
        cmin = 0
        cmax = 360
        norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)

        # Create a fake scalar mappable for the colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])  # This is necessary because we won't display the mappable

        # Create a new axis for the color wheel
        rect = [0.9, 0.1, 0.08, 0.5]
        ax_colorwheel = fig6.add_axes(rect, projection='polar', frame_on=False)

        # Create the color wheel data
        theta = np.linspace(0, 2 * np.pi, 512)  # Modified to fill the full circle
        radius = np.linspace(0.5, 1, 2)
        Theta, Radius = np.meshgrid(theta, radius)

        # Repeat the values from 0 to pi for the complete circle
        values = (Theta % np.pi) / np.pi  # Use modulo to repeat values

        ax_colorwheel.pcolormesh(Theta, Radius, values, cmap='hsv', shading='auto')

        # Make the color wheel visually appealing
        ax_colorwheel.set_yticks([])
        ax_colorwheel.set_theta_zero_location('E')
        ax_colorwheel.set_theta_direction(-1)
        fig6.delaxes(ax6[5])

        plt.show()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    files=[]
    main = MainWindow()
    main.show()
    sys.exit(app.exec())