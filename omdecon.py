#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Omdecon: deconvolution on the OMERO server. Omdecon is an OMERO wrapper for AIDA - Adaptative Image Deconvolution Algorithm 
by Étienne Labrie-Dion, 2017, Douglas Hospital Research Center, Montreal
"""
import time, os, glob, sys, shutil, tempfile

import omero, omero.scripts as scripts
from omero.gateway import BlitzGateway
from omero.rtypes import *

from PIL import Image
import numpy as np
import astropy.io.fits as fits
from astropy.io.fits.hdu.image import PrimaryHDU


#AIDA 1.41. For now Omdecon is not using the 3D capabilities of AIDA as 3D tif reading is broken in v1.41.
import aida.AIDA as AIDA
	
#Christoph Gohlke's PSF.py module for this (http://www.lfd.uci.edu/~gohlke/code/psf.py.html)
#Instructions to build the extension module available in psf.c
import psf

def process(conn, scriptParams):
	#The main part of the script, looping over all images in a Dataset or in the IDs list 
	
	
	#Read the script parameters
	datatype = scriptParams["Data_Type"]
	ids = scriptParams["IDs"]
	
	#Connect to OMERO and find the objects
	obs = conn.getObjects(datatype, ids)
	objects = list(obs)
	images = []
	if len(objects) == 0:
		print "Error: No %ss found for specified IDs" % datatype
		return
	if datatype == 'Dataset':
		for ds in objects:
			print "Processing Images from Dataset: %s" % ds.getName()
			imgs = list(ds.listChildren())
			images.extend(imgs)
	else:
		print "Processing Images identified by ID"
		images = objects	
	
	
	#Processing loop
	try: #The try/finally temporary directory code was taken from http://stackoverflow.com/questions/6884991/how-to-delete-dir-created-by-python-tempfile-mkdtemp on 2017/02/22
		
		imagepath = tempfile.mkdtemp() + "/" #create temp directory
		deconvoluted_images = [] #prepare the output list
		
		#Download all images
		for i in images:
			imageId = i.getId()
			sizeC = i.getSizeC()
			download(i, imageId, imagepath)
			get_PSF(i, imagepath, conn, scriptParams)
		
		#Deconvolve all images
		run_AIDA(images, imagepath) #AIDA is outside the loop because running AIDA more than once didn't work. AIDA supports dataset deconvolution so we just need to supply the list of images in run_AIDA().
		
		#Upload all images
		for i in images:
			newimage = upload_to_OMERO(i, imagepath, conn)
			deconvoluted_images.append(newimage) #Prepare the output
	
	finally: #remove temporary directory
		try:
			shutil.rmtree(imagepath)
		except OSError as exerror:
			if exerror.errno != errno.ENOENT:
				raise ("Couldn't remove the temporary directory. Please check the /tmp directory for leftover files.")
	
	#Format the output and return the images and messages to runAsScript
	if len(deconvoluted_images) == 0:
		return None, "No Images Deconvolved"
	if len(deconvoluted_images) == 1:
		msg = "Successfully Deconvoluted Image! New Image Name: %s" % deconvoluted_images[0].getName()
		return deconvoluted_images[0]._obj, msg
	else:
		ds = deconvoluted_images[0].getParent()
		if ds is not None:
			return ds._obj, "Successfully Deconvoluted %s Images from Dataset: %s!" % (len(deconvoluted_images),ds.getName())
		else:
			return None, "Successfully Deconvoluted %s Images!" % len(deconvoluted_images)

def get_PSF(stack, imagepath, conn, scriptParams):
	#Prepare the PSF for the deconvolution
	imageId = stack.getId()
	
	if scriptParams["Type_of_PSF"] == 'Measured PSF':
		#If a PSF is supplied, find the measured PSF and download it in the deconvolution folder.
		
		psfId = scriptParams["Measured_PSF_Image"]
		try:
			obs = conn.getObject("Image", psfId)
		except AttributeError:
			raise("The image supplied for the PSF (image %s) is not valid." %psfId)
		download(obs,str(imageId)+"PSF", imagepath)
	
	
	else:
		#If no measured PSF is supplied, the theoretical PSF must be calculated. 
		
		#Option1: Manually enter PSF parameters
		#Option2: Use OMERO metadata
		
		#First, prepare all the parameters that will be shared between the two options
		args = {}
		channelExcitation = []
		channelEmission = []
		sizeC = 0
		sizeX = stack.getPrimaryPixels().getSizeX()
		pixelX = stack.getPrimaryPixels().getPhysicalSizeX().getValue()
		
		#The PSF parameters are passed to psf.py through the args dictionary
		args["shape"] = (int(round(sizeX / 2)+1), int(round(sizeX / 2))+1) #Need to divide by 2 and add 1 to get the right final number of pixels
		args["dims"] = (int(round(pixelX * sizeX / 2)), int(round(pixelX * sizeX / 2))) #Because of the previous operation the PSF will be inconsequentially too small (1/x) and off-centered.
		args["pinhole_radius"] = 0.55 #default value
		args["pinhole_shape"] = 'round' #the Olympus FV1200 confocal's pinhole is square
		
		
		#Option 1: Manual PSF
		if scriptParams["Type_of_PSF"] == 'Manual theoretical PSF':
			
			fluoType = scriptParams["Fluorescence_Microscopy_"]
			
			args["pinhole_radius"] = scriptParams["Confocal_Pinhole______"]
			args["num_aperture"] = scriptParams["Objective_NA_____________"]
			
			if scriptParams["Imaging_Medium_________"] == 'Oil':
				args["refr_index"] = 1.515
			elif scriptParams["Imaging_Medium_________"] == 'Water':
				args["refr_index"] = 1.333
			elif scriptParams["Imaging_Medium_________"] == 'Air':
				args["refr_index"] = 1.0
			elif scriptParams["Imaging_Medium_________"] == 'CLARITY':
				args["refr_index"] = 1.45
				
			if scriptParams["_____Channel_1_-_Excitation"] !=0:
				sizeC += 1
				channelExcitation.append(scriptParams["_____Channel_1_-_Excitation"])
				channelEmission.append(scriptParams["_____Channel_1_-_Emission__"])
			
			if scriptParams["_____Channel_2_-_Excitation"] !=0:
				sizeC += 1
				channelExcitation.append(scriptParams["_____Channel_2_-_Excitation"])
				channelEmission.append(scriptParams["_____Channel_2_-_Emission__"])
			
			if scriptParams["_____Channel_3_-_Excitation"] !=0:
				sizeC += 1
				channelExcitation.append(scriptParams["_____Channel_3_-_Excitation"])
				channelEmission.append(scriptParams["_____Channel_3_-_Emission__"])
			
			if sizeC != stack.getSizeC():
				raise ValueError("Wrong number of PSF channels parameters set. Make sure to set the right number of channels.")
		
		
		#Option 2: Automatic PSF. (Uses the metadata from the stack. Whether this will work depends on the metadata.)
		else:	
		
			try:
				objective = stack.getObjectiveSettings().getObjective().getModel()
			except AttributeError:
				raise ("Couldn't identify the objective used for image %s. Please use the manual PSF settings." %imageId)

			try:
				args["num_aperture"] = stack.getObjectiveSettings().getObjective().getLensNA()
			except AttributeError:
				raise ("Couldn't find the numerical aperture for the objective used for image %s. Please use the manual PSF settings." %imageId)
				
			#Deduce immersion and microscopy type through the objective used:
			if objective in {"XLSLPLN25XGMP","XLPLN10XSVMP"}: #Two-photon CLARITY objectives
				fluoType = "TWOPHOTON"
				args["refr_index"] = 1.48
			elif objective in {"XLPLN25XWMP2", "XLUMPLFLN-W"}: #Two-photon water objectives #need to double check the metadata name of the XLUM objective
				fluoType = "TWOPHOTON"
				args["refr_index"] = 1.333
			elif objective == "PLAPON    60X O  NA:1.42": #Confocal oil objectives
				fluoType = "CONFOCAL"
				args["refr_index"] = 1.515
				args["pinhole_radius"] = 1 #TODO Figure out which metadata setting to select for pinhole_radius
			
				#TODO add air and widefield objectives	
			
			else:
				raise AttributeError("Couldn't identify the objective used for image %s. Please use the manual PSF settings." %imageId)
			
			for c in stack.getChannels():
				channel = c.getLogicalChannel()
				try:
					channelExcitation.append(channel.getExcitationWave().getValue())
					channelEmission.append(channel.getEmissionWave().getValue())
				except AttributeError:
					raise ("Couldn't find the excitation or emission wavelength of Channel %s from image %s. Please use the manual PSF settings." %(c,imageID))
			
			sizeC = stack.getSizeC()			
		
		#Calculate the PSF.
		for c in range(sizeC): #Loop through all channels since each one has a different theoretical PSF
		
			args["ex_wavelen"] = channelExcitation[c]
			args["em_wavelen"] = channelEmission[c]
			
			#This is the step where the image of the theoretical PSF is produced
			if fluoType == "TWOPHOTON":
				obsvol = psf.PSF(psf.GAUSSIAN | psf.TWOPHOTON, **args) #The settings are supplied through the args dictionary
			elif fluoType == "CONFOCAL":
				obsvol = psf.PSF(psf.GAUSSIAN | psf.CONFOCAL, **args)
			elif fluoType == "WIDEFIELD":
				obsvol = psf.PSF(psf.GAUSSIAN | psf.WIDEFIELD, **args)
			
			#The PSF needs to be completed with the mirror_symmetry function and corrected because it's produced with 
			#one extra row and one extra column of pixels; we need to shave off those two extra lines. This doesn't seem ideal but this 
			#is the recommandation from C. Gohlke and AIDA centers the PSF anyway.
			image_psf = psf.mirror_symmetry(obsvol.slice(0)) 
			image_psf = np.delete(image_psf,1,0) 
			image_psf = np.delete(image_psf,1,1) 
			
			#Donwload the PSF. This is the same operation as the download function
			#Convert to FITS and save
			pixels = image_psf.astype(np.float32)
			hdu = PrimaryHDU(pixels)
			hdulist = fits.HDUList([hdu])
			deconv_path = os.path.join(imagepath, "%s%s_C%s_%03d_%03d.fits" % (imageId,"PSF",c,0,0))
			hdulist.writeto(deconv_path)

def download(stack, imageId, imagepath):
	#Download the Images to the temp folder
	print " Downloading %s: %s" % (imageId,stack.name)
	
	sizeZ = stack.getSizeZ()
	sizeC = stack.getSizeC()
	sizeT = stack.getSizeT()	
	zctList = [] 
	
	if stack.getPrimaryPixels().getSizeX() != stack.getPrimaryPixels().getSizeY():
		raise("The image %s is not square. Please supply a square image only." %imageId)
			
	for z in range(sizeZ): #build the zctList as described in the OMERO documentation
		for c in range(sizeC):
			for t in range(sizeT):
				zctList.append((z, c, t))
	
	pixels = stack.getPrimaryPixels()
	planes = pixels.getPlanes(zctList) #planes is the OMERO generator (images are loaded only when they are accessed) for the images.
	
	for index, plane in enumerate(planes): #Iterate through all planes
		
		plane = plane.astype(np.float32) #AIDA needs 32bit images
		
		if "PSF" not in str(imageId) :
			#measure background if not a PSF image
			background = int(plane.min()+plane.std()*.25)
			if background == 0: #if the whole image is blank, background is 1
				background = 1
		else:
			background = 0
		
		plane = np.subtract(plane, background) #AIDA needs negative pixel values for the background
		
		#Convert to FITS and save
		hdu = PrimaryHDU(plane)
		hdulist = fits.HDUList([hdu])
		deconv_path = os.path.join(imagepath, "%s_C%s_%03d_%03d.fits" % (imageId,zctList[index][1],zctList[index][2],zctList[index][0]))
		hdulist.writeto(deconv_path)

def run_AIDA(images,imagepath):

	omeroSettings = imagepath + 'omeroSettings.txt'
	path = imagepath
	results_directory = path + 'Results'
	image_filenames = []
	PSF_filenames = []
	dataset_label = []
	decon_type = "'myopic'" 
	
	#Build the list of images to be deconvoluted
	for i in images:
		imageId = i.getId()
		c = i.getSizeC()
		
		for channel in range(c):
			image_filenames.append('%s%s_C%s_000_' % (imagepath,imageId,channel))
			PSF_filenames.append('%s%sPSF_C%s_000_' % (imagepath,imageId,channel))
			dataset_label.append('%s_C%s' % (imageId,channel))
			
		
	with open(omeroSettings, 'wb') as settings_file:
		settings_file.write("path = '%s' \nresults_directory = '%s' \nimage_filenames = %s \nPSF_filenames = %s \ndataset_label = %s \ndecon_type = %s" % (path,results_directory,image_filenames,PSF_filenames,dataset_label,decon_type))
	
	print "Start of Deconvolution"
	AIDA.RunAIDA(omeroSettings)

	#Restore the stdout from the AIDA log, otherwise all of the print statements are caught by the AIDA log
	sys.stdout = sys.__stdout__	

def upload_to_OMERO(stack, imagepath, conn):
	print "Preparing to upload the deconvolution of image '%s'" % (stack.name)
	
	imageName = stack.name
	datasetID = stack.getParent()
	imageId = stack.getId()
	sizeZ = stack.getSizeZ()
	sizeC = stack.getSizeC()
	sizeT = stack.getSizeT()
	pixelType = stack.getPixelsType()
	clist = range(sizeC)
	zctList = [] 
	newName = imageName[:-4] + "_Deconv" + imageName[-4:]
	
	for z in range(sizeZ): 
		for c in range(sizeC):
			for t in range(sizeT):
				zctList.append((z, c, t))

	deconvoluted_path = glob.glob(imagepath + 'Results/' + "%s_C*" % (imageId))
	print "Found deconvoluted images folder starting at '%s'" % (deconvoluted_path[0])
	
	def planeGen():
		for i in zctList:
			find_image = glob.glob(deconvoluted_path[i[1]] + "/Robj*%s_C%s_%03d_%03d*" % (imageId,i[1],i[2],i[0]))[0]
			open_image = fits.open(find_image)
			plane = open_image[0].data
			plane = plane.astype(pixelType)
			yield plane
	
	print "Uploading deconvoluted Image '%s'" % newName
	
	#uploads the deconvoluted image to OMERO and automatically copies the metadata by supplying the sourceImageId	
	uploaded = conn.createImageFromNumpySeq(planeGen(), newName, sizeZ, sizeC, sizeT, description="Image Deconvoluted with AIDA", dataset = datasetID, sourceImageId = imageId, channelList = clist)
	
	print "Created new image, ID: %s, Name: %s" % (uploaded.getId(),newName)
	
	return uploaded
	
def runAsScript():

	dataTypes = [rstring('Dataset'), rstring('Image')]
	
	microscopyTypes = [rstring('CONFOCAL'), rstring('TWOPHOTON'), rstring('WIDEFIELD')]
	
	immersionTypes = [rstring('Oil'), rstring('Water'), rstring('Air'), rstring('CLARITY')]
	
	psfTypes = [rstring('Measured PSF'), rstring('Manual theoretical PSF'), rstring('Automatic theoretical PSF')]
	
	client = scripts.client('omdecon',
	"""omdecon: deconvolution on the OMERO server.
	
	Omdecon uses AIDA, the Adaptative Image Deconvolution Algorithm, to deconvolve images. 
	More info: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3166524/
	Also: https://github.com/erikhom/aida ; https://code.google.com/archive/p/aida-deconvolution/
	
	
	Omdecon can calculate the theoretical PSF if you don't have one on hand. This feature uses Christoph Gohlke's PSF.py module.
	More info: http://www.lfd.uci.edu/~gohlke/code/psf.py.html
	
	READ BEFORE YOU USE OMDECON:
		1- Omdecon only works with square images (i.e. 512x512, 600x600).
		2- The automatic theoretical PSF calculation requires the following metadata:
		   -Objective, Numerical Aperture, Excitation and Emission for all channels and Pixel Calibration.
		3- You need to input the right number of channels when setting the PSF manually.
		4- Images are saved in the same dataset with "_Deconv" added to the name.
	
	""",
		
		scripts.String("Data_Type", optional=False, grouping="1", 
		description="Type of data to deconvolve", values=dataTypes, default=rstring('Image')),
		
		scripts.List("IDs", optional=False, grouping="1.1",
		description = "List of images to deconvolve.").ofType(rlong(0)),
				
		scripts.String("PSF Configuration", grouping = "2",
		description="PSF settings", values=None, default=None),
		
		scripts.String("Type_of_PSF", grouping = "2.1",
		description="Choose the type of PSF you want to use.", values=psfTypes, default=rstring('Measured PSF')),
		
		scripts.Long("Measured_PSF_Image", optional=True, grouping="2.2",
		description="Input '0' to calculate the theoretical PSF", default=0),
		
		scripts.String("Manual settings for theoretical PSF (ignored for automatic)", grouping = "3",
		description="Supply the parameters here to generate a theoretical PSF with manual settings", values=None, default=None),
		
		scripts.String("Fluorescence_Microscopy_", grouping="3.01",
		description="Choose the type of fluorescence microscopy used", values=microscopyTypes, default=rstring('CONFOCAL')),
		
		scripts.String("Imaging_Medium_________", grouping="3.02",
		description="Select the media used for imaging", values=immersionTypes, default=rstring('Oil')),
		
		scripts.Float("Objective_NA_____________", optional=False, grouping="3.03",
		description="Supply the numerical aperture of the objective", default=1.42),
		
		scripts.Float("Pixel_Size_in_microns", optional=False, grouping="3.04",
		description="Enter the size of your pixels", default=0.05),
		
		scripts.Float("Confocal_Pinhole______", optional=False, grouping="3.05",
		description="Enter the pinhole radius for confocal microscopy", default=0.55),
		
		scripts.Int("_____Channel_1_-_Excitation",default=488, grouping="4.1", 
		description="Input '0' if not using this channel"),
		
		scripts.Int("_____Channel_1_-_Emission__",default=510, grouping="4.2", 
		description="Input '0' if not using this channel"),
			
		scripts.Int("_____Channel_2_-_Excitation",default=0, grouping="4.3", 
		description="Input '0' if not using this channel"),
		
		scripts.Int("_____Channel_2_-_Emission__",default=0, grouping="4.4", 
		description="Input '0' if not using this channel"),		

		scripts.Int("_____Channel_3_-_Excitation",default=0, grouping="4.5", 
		description="Input '0' if not using this channel"),
		
		scripts.Int("_____Channel_3_-_Emission__",default=0, grouping="4.6", 
		description="Input '0' if not using this channel"),
		
		
		version="1.1.0",
		authors=["Etienne Labrie-Dion"],
		institutions=["Douglas Hospital Research Center, McGill University, Montreal"],
		contact="etienne.labrie-dion@douglas.mcgill.ca",		
	)
	
	try:
		#Obtain the user parameters
		scriptParams = {}
		conn = BlitzGateway(client_obj=client)

		client.enableKeepAlive(240)
		
		for key in client.getInputKeys():
				if client.getInput(key):
					scriptParams[key] = client.getInput(key, unwrap=True)
		print scriptParams
		
		robj, message = process(conn, scriptParams)
		client.setOutput("Message", rstring(message))
		
		if robj is not None:
				client.setOutput("Result", robject(robj))

	finally:
	
	#TODO: send an email confirming the end of the deconvolution and success/failure
		client.closeSession()

if __name__ == "__main__":
	runAsScript()