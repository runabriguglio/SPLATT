"""
Authors
  - C. Selmi:  written in 2019
               modified in 2021
"""

import os
import logging
from astropy.io import fits as pyfits
import numpy as np
from m4.configuration import config_folder_names as fold_name
from m4.configuration.ott_parameters import OttParameters
from m4.ground import tracking_number_folder
from m4.ground import zernike
from m4.configuration.ott_parameters import OtherParameters

WHO_PAR_AND_RM = "PAR + RM"
WHO_PARABLE = "PAR"
WHO_RM = "RM"
WHO_M4 = "M4"


class OpticalCalibration:
    """
    Class for the optical calibration and interaction matrix creation

    HOW TO USE IT::

        from m4.utils.optical_calibration import OpticalCalibration
        cal = OpticalCalibration(ott, interf)
        cal.measureAndAnalysisCalibrationMatrix(who, command_amp_vector,
                                            n_push_pull, n_frames, delay)

    """

    def __init__(self, ott, interf):
        """The constructor"""
        self._logger = logging.getLogger("OPT_CALIB:")
        self._interf = interf
        self._ott = ott
        # start
        self._nPushPull = None
        self._commandAmpVector = None
        self._who = None
        # from calibration
        self._dofIndex = None
        self._commandMatrix = None
        self._commandList = None
        self.tt = None
        # from analysis
        self._cube = None
        self._mask = None
        self._intMat = None

        self._fullCommandMatrix = None
        self._fullCube = None

    @staticmethod
    def _storageFolder():
        """Creates the path where to save measurement data"""
        return fold_name.CALIBRATION_ROOT_FOLDER

    def measureAndAnalysisCalibrationMatrix(
        self, who, command_amp_vector, n_push_pull, n_frames, delay, tnpar
    ):
        """
        Parameters
        ----------
        who: string
            string indicating the optical element
            on which to perform the calibration
            cal.WHO_PAR_AND_RM for parabola and reference mirror
            cal.WHO_PARABLE for parabola (non implemented)
            cal.WHO_RM for reference mirror (not implemented)
            cal.WHO_M4 for deformable mirror
        command_amp_vector: numpy array [mm]
                            vector containing the amplitude of the
                            commands to give degrees of freedom to
                            calibrate
        n_push_pull: int
                    number of push pull
        n_frames: int
                number of frame for 4D measurement

        Returns
        -------
        tt : string
            tracking number containing measurements and IntMat from analysis
        """
        self._nPushPull = n_push_pull
        self._commandAmpVector = command_amp_vector

        dove, self.tt = tracking_number_folder.createFolderToStoreMeasurements(
            self._storageFolder()
        )
        self._logger.info("Measure of calibration. Location: %s", self.tt)

        # measurement
        dofIndex_vector = self._logAndDefineDovIndexForCommandMatrixCreation(
            who)
        self._commandMatrix, self._commandList = self.createCmatAndCmdList(
            command_amp_vector, dofIndex_vector
        )
        self._measureAndStore(self._commandList, dove, n_frames, delay)
        # analysis
        self._createCube(False)  # cubo non normalizzato
        cube = self.getCube()
        self._mask = self._findMask(cube)
        self._intMat = self.getInteractionMatrix(tnpar)

        self._saveCalibrationInfoAndResultsAsFits(dove)

        #         if self._mixed_method == False:
        #             self._createCube() #cubo normalizzato
        #             cube_norm = self.getCube()
        #             self._mask = self._findMask(cube_norm)
        #             self._intMatNorm = self.getInteractionMatrix(self._mask)
        #             fits_file_name = os.path.join(dove,
        # 'InteractionMatrixNorm.fits')
        #             pyfits.writeto(fits_file_name, self._intMatNorm,
        # overwrite=True)
        return self.tt

    def _findMask(self, cube):
        ima = cube[:, :, 0]
        from m4.utils import roi

        rois = roi.roiGenerator(ima)
        if fold_name.simulated_interf is True:
            mask_index = OtherParameters.MASK_INDEX_SIMULATORE
        else:
            mask_index = OtherParameters.MASK_INDEX_TOWER
        mask = rois[mask_index]
        return mask

    def getCommandMatrix(self):
        """
        Returns
        -------
        commandMatrix: numpy array
                    command matrix used for calibration
        """
        return self._commandMatrix

    def getFullCommandMatrix(self):
        """
        Returns
        -------
        commandMatrix: numpy array
                    command matrix used for calibration
        """

        # Split each array into a list of columns
        arr1 = self._commandMatrix.copy()
        arr2 = -arr1
        columns1 = np.hsplit(arr1, arr1.shape[1])
        columns2 = np.hsplit(arr2, arr2.shape[1])

        # Alternate columns from each array
        interlaced_columns = [
            column for pair in zip(columns1, columns2) for column in pair
        ]

        # Horizontally stack the alternated columns back into a single array
        interlaced_array = np.hstack(interlaced_columns)

        self._fullCommandMatrix = np.tile(interlaced_array, self._nPushPull)

        return self._fullCommandMatrix

    def getMask(self):
        """
        Returns
        -------
        mask: numpy array
            mask used for interaction matrix calculation
        """
        return self._mask

    def getWho(self):
        return self._who

    def _measureAndStore(self, command_list, dove, n_frames, delay):
        if self._who == "PAR + RM":
            vec_push_pull = np.array((1, -1))
            # mis = (len(command_list)-2) * n_push_pull * vec_push_pull.size
            par0 = self._ott.parabola.getPosition()
            rm0 = self._ott.referenceMirror.getPosition()
            for k in range(self._nPushPull):
                for i in range(len(command_list) - 2):
                    j = (len(command_list) - 2) * k * 2
                    mis = np.array([j, j + 1])
                    if i == 0:
                        pcmd = np.array(command_list[i])
                        for v in range(vec_push_pull.size):
                            par1 = pcmd * vec_push_pull[v]
                            print(par1)
                            self._ott.parabola.setPosition(par0 + par1)
                            masked_ima = self._interf.acquire_phasemap(
                                n_frames, delay)
                            masked_ima = self._interf.intoFullFrame(masked_ima)

                            name = "Frame_%04d.fits" % (2 * i + mis[v])
                            print(name)
                            self._interf.save_phasemap(dove, name, masked_ima)
                            self._ott.parabola.setPosition(par0)
                    elif i == 1 or i == 2:
                        if i == 1:
                            l = i
                        else:
                            l = i + 1
                        pcmd = np.array(command_list[l])
                        rcmd = np.array(command_list[l + 1])
                        for v in range(vec_push_pull.size):
                            par1 = pcmd * vec_push_pull[v]
                            rm1 = rcmd * vec_push_pull[v]
                            self._ott.parabola.setPosition(par0 + par1)
                            if np.count_nonzero(rm1) != 0:
                                self._ott.referenceMirror.setPosition(
                                    rm0 + rm1)
                            print(par1, rm1)
                            masked_ima = self._interf.acquire_phasemap(
                                n_frames, delay)
                            masked_ima = self._interf.intoFullFrame(masked_ima)

                            name = "Frame_%04d.fits" % (2 * i + mis[v])
                            print(name)
                            self._interf.save_phasemap(dove, name, masked_ima)
                            self._ott.parabola.setPosition(par0)
                            self._ott.referenceMirror.setPosition(rm0)
                    else:
                        rcmd = np.array(command_list[i + 2])
                        for v in range(vec_push_pull.size):
                            rm1 = rcmd * vec_push_pull[v]
                            self._ott.referenceMirror.setPosition(rm0 + rm1)
                            print(rm1)
                            masked_ima = self._interf.acquire_phasemap(
                                n_frames, delay)
                            masked_ima = self._interf.intoFullFrame(masked_ima)
                            name = "Frame_%04d.fits" % (2 * i + mis[v])
                            print(name)
                            self._interf.save_phasemap(dove, name, masked_ima)
                            self._ott.referenceMirror.setPosition(rm0)
        elif self._who == "PAR":
            pass
        elif self._who == "RM":
            pass
        elif self._who == "M4":
            vec_push_pull = np.array((1, -1))
            m0 = self._ott.m4Exapode.getPosition()
            for k in range(self._nPushPull):
                for i in range(len(command_list)):
                    j = (len(command_list)) * k * 2
                    mis = np.array([j, j + 1])
                    for v in range(vec_push_pull.size):
                        m4_cmd = command_list[i] * vec_push_pull[v]
                        self._ott.me.setPosition(m0 + m4_cmd)
                        masked_ima = self._interf.acquire_phasemap(
                            n_frames, delay)
                        masked_ima = self._interf.intoFullFrame(masked_ima)
                        name = "Frame_%04d.fits" % (2 * i + mis[v])
                        self._interf.save_phasemap(dove, name, masked_ima)
            self._ott.m4Exapode.setPosition(m0)

    def _logAndDefineDovIndexForCommandMatrixCreation(self, who):
        """
        who:
            cal.WHO_PAR_AND_RM for parabola and reference mirror
            cal.WHO_PARABLE for parabola (non implemented)
            cal.WHO_RM for reference mirror (not implemented)
            cal.WHO_M4 for deformable mirror
        """
        if who == "PAR + RM":
            self._who = "PAR + RM"
            self._dofIndex = np.append(
                OttParameters.PARABOLA_DOF, OttParameters.RM_DOF)
        elif who == "PAR":
            self._who = "PAR"
            self._dofIndex = OttParameters.PARABOLA_DOF
        elif who == "RM":
            self._who = "RM"
            self._dofIndex = OttParameters.RM_DOF
        elif who == "M4":
            self._who = "M4"
            self._dofIndex = OttParameters.M4_DOF
        else:
            raise OSError("Who= %s doesnt exists" % who)

        self._logger.info("Creation of the command matrix for %s", self._who)
        return self._dofIndex

    def createCmatAndCmdList(self, command_amp_vector, dofIndex_vector):
        """
        Function to allow the creation of the matrix of commands and the
        decomposition of them in a list of commands to assign to devices

        Parameters
        ----------
        command_amp_vector: numpy array [mm]
            vector containing the amplitude of the
            commands to give degrees of freedom to calibrate
        dofIndex_vector: numpy array
            vector containing position of Dof to move in standard vector of
            six position for command devices

        Returns
        -------
        command_matrix: numpy array
            matrix 5x5 composed using command_amp_vector values and relationship between them
        command_list: list
            decomposition of command matrix in list of command
        """
        # crea matrice 5 x 5
        command_matrix = np.zeros(
            (command_amp_vector.size, command_amp_vector.size))
        for i in range(command_amp_vector.shape[0]):
            j = i
            if i == 1 or i == 2:
                command_matrix[i, j] = command_amp_vector[i]
                command_matrix[i, j + 2] = (
                    OttParameters.par_rm_coef_for_coma_measuremets
                    * command_amp_vector[i]
                )
            else:
                command_matrix[i, j] = command_amp_vector[i]
        #         elif mixed_method == False:
        #             command_matrix = np.zeros((command_amp_vector.size,
        #  command_amp_vector.size))
        #             for i in range(command_amp_vector.shape[0]):
        #                 command_matrix[i, i] = command_amp_vector[i]

        # crea i comandi
        command_list = []
        for i in range(command_matrix.shape[0]):
            if i == 1 or i == 2:
                cmd = np.zeros(6)
                cmd[dofIndex_vector[i]] = command_matrix[i, i]
                command_list.append(cmd)
                cmd1 = np.zeros(6)
                cmd1[dofIndex_vector[i + 2]] = command_matrix[i, i + 2]
                command_list.append(cmd1)
            else:
                # cmd_amp = command_matrix[i, i]
                cmd = np.zeros(6)
                cmd[dofIndex_vector[i]] = command_matrix[i, i]
                command_list.append(cmd)
        return command_matrix, command_list

    def _saveCalibrationInfoAndResultsAsFits(self, dove):
        """
        Save fits file for the command matrix and the data relating to
        its creation

        args:
            dove = path that indicates where to save the command matrix file
        """
        fits_file_name = os.path.join(dove, "CalibrationInfo.fits")
        header = pyfits.Header()
        header["NPUSHPUL"] = self._nPushPull
        header["WHO"] = self._who
        pyfits.writeto(fits_file_name, self._commandAmpVector, header)
        pyfits.append(fits_file_name, self._commandMatrix.T, header)
        pyfits.append(fits_file_name, self._mask.astype(int), header)
        pyfits.append(fits_file_name, self._intMat, header)

        # file separti per Runa
        fits_file_name = os.path.join(dove, "CommandAmplitude.fits")
        pyfits.writeto(fits_file_name, self._commandAmpVector)
        fits_file_name = os.path.join(dove, "CMat.fits")
        pyfits.writeto(fits_file_name, self._commandMatrix.T)
        fits_file_name = os.path.join(dove, "Mask.fits")
        pyfits.writeto(fits_file_name, self._mask.astype(int))
        fits_file_name = os.path.join(dove, "InteractionMatrix.fits")
        pyfits.writeto(fits_file_name, self._intMat)

    @staticmethod
    def loadCalibrationObjectFromFits(tt):
        """Creates the object using information contained in calibration fits file

        Parameters
        ----------
        tt: string
            tracking number

        Returns
        -------
        theObject: object
                 opt_calibration class object
        """
        ott = None
        interf = None
        theObject = OpticalCalibration(ott, interf)
        theObject.tt = tt
        dove = os.path.join(theObject._storageFolder(), tt)
        file = os.path.join(dove, "CalibrationInfo.fits")
        header = pyfits.getheader(file)
        hduList = pyfits.open(file)
        theObject._who = header["WHO"]
        theObject._nPushPull = header["NPUSHPUL"]
        theObject._commandAmpVector = hduList[0].data
        theObject._commandMatrix = hduList[1].data
        theObject._mask = hduList[2].data
        theObject._intMat = hduList[3].data
        return theObject

    def _createFullCube(self, norm=True):
        """ """
        if norm != True:
            self._commandAmpVector = np.ones(self._commandAmpVector.size)
        self._logger.info("Creation of the cube relative to %s", self.tt)
        self._fullCube = None
        self._fold = os.path.join(OpticalCalibration._storageFolder(), self.tt)

        if self._fullCommandMatrix is None:
            dummy = self.getFullCommnadMatrix()
        for i in range(self._fullCommandMatrix.shape[1]):
            name_pos = "Frame_%04d.fits" % i
            print("Reding " + name_pos)
            file = os.path.join(self._fold, name_pos)
            hduList = pyfits.open(file)
            final_ima = np.ma.masked_array(
                hduList[0].data, mask=hduList[1].data.astype(bool)
            )

            if self._fullCube is None:
                self._fullCube = final_ima
            else:
                self._fullCube = np.ma.dstack((self._fullCube, final_ima))
        return

    def _createCube(self, norm=True):
        """ """
        if norm != True:
            self._commandAmpVector = np.ones(self._commandAmpVector.size)
        self._logger.info("Creation of the cube relative to %s", self.tt)
        self._cube = None
        self._fold = os.path.join(OpticalCalibration._storageFolder(), self.tt)
        for i in range(self._commandAmpVector.shape[0]):
            for j in range(self._nPushPull):
                k = 2 * i + 2 * self._commandAmpVector.shape[0] * j
                name_pos = "Frame_%04d.fits" % k
                name_neg = "Frame_%04d.fits" % (k + 1)
                file = os.path.join(self._fold, name_pos)
                hduList = pyfits.open(file)
                image_pos = np.ma.masked_array(
                    hduList[0].data, mask=hduList[1].data.astype(bool)
                )
                file = os.path.join(self._fold, name_neg)
                hduList = pyfits.open(file)
                image_neg = np.ma.masked_array(
                    hduList[0].data, mask=hduList[1].data.astype(bool)
                )

                if self._who == "PAR + RM":
                    image = (image_pos - image_neg) / \
                        (2 * self._commandAmpVector[i])
                elif self._who == "M4":
                    image = (image_pos + image_neg) / \
                        (2 * self._commandAmpVector[i])

                if j == 0:
                    all_push_pull_act = image
                else:
                    all_push_pull_act = np.ma.dstack(
                        (all_push_pull_act, image))
            if self._nPushPull == 1:
                final_ima = all_push_pull_act
            else:
                final_ima = np.ma.mean(all_push_pull_act, axis=2)

            if self._cube is None:
                self._cube = final_ima
            else:
                self._cube = np.ma.dstack((self._cube, final_ima))
        return

    def getCube(self):
        """
        Returns
        -------
        cube: numpy masked array
            analyzed measurements
        """
        if self._cube is None:
            self._createCube(False)  # lo crea sempre non normalizzato
            self._cube = self.getCube()
        return self._cube

    def getFullCube(self):
        """
        Returns
        -------
        cube: numpy masked array
            analyzed measurements
        """
        if self._fullCube is None:
            self._createFullCube(False)  # lo crea sempre non normalizzato
        return self._fullCube

    def _createInteractionMatrix(self, mask, tnpar):
        coefList = []
        self._cube = self.getCube()
        for i in range(self._cube.shape[2]):
            ima = np.ma.masked_array(self._cube[:, :, i], mask=mask)
            coef, mat = zernike.zernikeFit(ima, np.arange(10) + 1)

            # modRB20231026 to implement aux mask fitting. the following lines replace the previous one
            # from m4.utils import image_registration_lib as imgreg
            from m4.ground import geo

            # img = self._interf.intoFullFrame(ima) this was removed since the frames are already saved in fullframe
            ### now is passed! modRB 20240725 tnpar = "20240521_161525"  # '20231016_124531'
            if tnpar is not None:
                print("Using global modes fitting, TNPar: " + tnpar)
                # par = imgreg.load_registeredPar(tnpar)
                par = self._load_registeredPar(tnpar)

                cir = geo.qpupil(-1 * par.mask + 1)
                mm = geo.draw_mask(
                    par.data * 0, cir[0], cir[1], 1.44 / 0.00076 / 2, out=0)
                coef, mat = zernike.zernikeFitAuxmask(
                    ima, mm, np.arange(10) + 1)  # was img
            # end modRB
            else:
                coef, mat = zernike.zernikeFit(ima, np.arange(10) + 1)
            # z= np.array([2,3,4,7,8])
            z = np.array([1, 2, 3, 6, 7])
            final_coef = np.zeros(z.shape[0])
            final_coef = coef[z]
            self._mat = mat
            coefList.append(final_coef)

        self._intMat = np.zeros((coefList[0].shape[0], self._cube.shape[2]))
        for j in range(self._cube.shape[2]):
            self._intMat.T[j] = coefList[j]

    def getFullLocalInteractionMatrix(self):
        coefList = []
        self._cube = self.getFullCube()
        for i in range(self._fullCube.shape[2]):
            ima = np.ma.masked_array(self._fullCube[:, :, i], mask=self._mask)
            coef, mat = zernike.zernikeFit(ima, np.arange(10) + 1)

            # z= np.array([2,3,4,7,8])
            z = np.array([1, 2, 3, 6, 7])
            final_coef = np.zeros(z.shape[0])
            final_coef = coef[z]
            self._mat = mat
            coefList.append(final_coef)

        self._fullIntMat = np.zeros(
            (coefList[0].shape[0], self._cube.shape[2]))
        for j in range(self._fullCube.shape[2]):
            self._fullIntMat.T[j] = coefList[j]
        return self._fullIntMat

    def getInteractionMatrix(self, tnpar=None):
        """
        Returns
        -------
        intMat: numpy array
                interaction matrix
        """
        if self._intMat is None:
            self._createInteractionMatrix(self._mask, tnpar)
        return self._intMat

    def _load_registeredPar(self, tn, zlist=[1, 2, 3, 4]):
        # fold=th.foldname.PARABOLA_REMAPPED_FOLDER+'/'+tn+'/'
        name = fold_name.PARABOLA_REMAPPED_FOLDER + "/" + tn + "/" + "par_remapped.fits"
        print("Loading registered Par " + name)
        hdu = pyfits.open(name)
        img = hdu[0].data
        mask = hdu[1].data
        imgout = np.ma.masked_array(img, mask)
        coeff, mat = zernike.zernikeFit(imgout, zlist)
        surf = zernike.zernikeSurface(imgout, coeff, mat)
        imgout = imgout - surf
        return imgout
