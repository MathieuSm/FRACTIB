# 00 Script initialization
import os
from paraview.simple import *

## Disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()



# 01 Define functions
def PrintSlice(MHDImage,RenderView_1,SaveDir,ImageOpacity=1):

    # create a new 'Meta File Series Reader'
    Image_1 = MetaFileSeriesReader(registrationName=MHDImage[:-4], FileNames=[SaveDir+MHDImage])
    ImageDisplay = Show(Image_1, RenderView_1, 'UniformGridRepresentation')
    
    # Get central position for slicing
    ImageSize = Image_1.GetDataInformation().GetExtent()
    X_SlicePosition = ImageSize[1] / 2
    Y_SlicePosition = ImageSize[3] / 2
    Z_SlicePosition = ImageSize[5] / 2
    SlicePosition = [X_SlicePosition, Y_SlicePosition, Z_SlicePosition]

    # hide data in view
    Hide(Image_1, RenderView_1)

    # create a new 'Slice'
    Slice_1 = Slice(registrationName=MHDImage[:-4]+ '_Slice', Input=Image_1)
    Slice_1.SliceType = 'Plane'
    Slice_1.HyperTreeGridSlicer = 'Plane'
    Slice_1.SliceOffsetValues = [0.0]

    # init the 'Plane' selected for 'SliceType'
    Slice_1.SliceType.Origin = [X_SlicePosition, Y_SlicePosition, Z_SlicePosition]

    # init the 'Plane' selected for 'HyperTreeGridSlicer'
    Slice_1.HyperTreeGridSlicer.Origin = [X_SlicePosition, Y_SlicePosition, Z_SlicePosition]

    # show data in view
    Slice_1Display = Show(Slice_1, RenderView_1, 'GeometryRepresentation')

    # get color transfer function/color map for 'MetaImage'
    MetaImageLUT = GetColorTransferFunction('MetaImage')
    MetaImageLUT.RGBPoints = [-1699.5, 0.231373, 0.298039, 0.752941, 6326.9501953125, 0.865003, 0.865003, 0.865003, 14353.400390625, 0.705882, 0.0156863, 0.14902]
    MetaImageLUT.ColorSpace = 'RGB'
    MetaImageLUT.ScalarRangeInitialized = 1.0

    # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
    MetaImageLUT.ApplyPreset('X Ray', True)

    # get color legend/bar for MetaImageLUT in view RenderView_1
    MetaImageLUTColorBar = GetScalarBar(MetaImageLUT, RenderView_1)
    MetaImageLUTColorBar.AutoOrient = 0
    MetaImageLUTColorBar.Orientation = 'Horizontal'
    MetaImageLUTColorBar.WindowLocation = 'LowerCenter'
    MetaImageLUTColorBar.ComponentTitle = ''

    # Properties modified on MetaImageLUTColorBar
    MetaImageLUTColorBar.Title = 'Gray values (-)'
    MetaImageLUTColorBar.TitleFontFamily = 'Times'
    MetaImageLUTColorBar.TitleFontSize = 24
    MetaImageLUTColorBar.LabelFontFamily = 'Times'
    MetaImageLUTColorBar.LabelFontSize = 24

    # Properties modified on MetaImageLUTColorBar
    MetaImageLUTColorBar.AddRangeLabels = 0

    # trace defaults for the display properties.
    Slice_1Display.Representation = 'Surface'
    Slice_1Display.ColorArrayName = ['POINTS', 'MetaImage']
    Slice_1Display.LookupTable = MetaImageLUT
    Slice_1Display.SelectTCoordArray = 'None'
    Slice_1Display.SelectNormalArray = 'None'
    Slice_1Display.SelectTangentArray = 'None'
    Slice_1Display.OSPRayScaleArray = 'MetaImage'
    Slice_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    Slice_1Display.SelectOrientationVectors = 'None'
    Slice_1Display.ScaleFactor = 47.400000000000006
    Slice_1Display.SelectScaleArray = 'MetaImage'
    Slice_1Display.GlyphType = 'Arrow'
    Slice_1Display.GlyphTableIndexArray = 'MetaImage'
    Slice_1Display.GaussianRadius = 2.37
    Slice_1Display.SetScaleArray = ['POINTS', 'MetaImage']
    Slice_1Display.ScaleTransferFunction = 'PiecewiseFunction'
    Slice_1Display.OpacityArray = ['POINTS', 'MetaImage']
    Slice_1Display.OpacityTransferFunction = 'PiecewiseFunction'
    Slice_1Display.DataAxesGrid = 'GridAxesRepresentation'
    Slice_1Display.PolarAxes = 'PolarAxesRepresentation'
    
    # Modify opacity
    Slice_1Display.Opacity = ImageOpacity

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    Slice_1Display.OSPRayScaleFunction.Points = [-82.39119720458984, 0.0, 0.5, 0.0, 82.35810089111328, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    Slice_1Display.ScaleTransferFunction.Points = [-1699.5, 0.0, 0.5, 0.0, 14097.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    Slice_1Display.OpacityTransferFunction.Points = [-1699.5, 0.0, 0.5, 0.0, 14097.0, 1.0, 0.5, 0.0]

    # show color bar/color legend
    Slice_1Display.SetScalarBarVisibility(RenderView_1, True)

    # set active source
    SetActiveSource(Image_1)

    # update the view to ensure updated data information
    RenderView_1.Update()


    # get layout
    layout1 = GetLayout()

    # layout/tab size in pixels
    layout1.SetSize(998, 522)

    # save screenshot
    SaveScreenshot(SaveDir+'00_'+MHDImage[:-4]+'_YZPlane.png', RenderView_1, ImageResolution=[998, 522])

    return Image_1, Slice_1, Slice_1Display, SlicePosition

def PrintDeformations(SliceDisplay,SlicePosition,ScalarValues_1,Min,Max,RenderView_1,SaveDir,ImageOpacity=0.5):
    
    # create a new 'Meta File Series Reader'
    Image_1 = MetaFileSeriesReader(registrationName=ScalarValues_1[:-4], FileNames=[SaveDir+ScalarValues_1])
    
    # Get central position for slicing
    X_SlicePosition = SlicePosition[0] - 0.01
    Y_SlicePosition = SlicePosition[1]
    Z_SlicePosition = SlicePosition[2]

    # hide data in view
    Hide(Image_1, RenderView_1)

    # create a new 'Slice'
    Slice_1 = Slice(registrationName=ScalarValues_1[:-4]+ '_Slice', Input=Image_1)
    Slice_1.SliceType = 'Plane'
    Slice_1.HyperTreeGridSlicer = 'Plane'
    Slice_1.SliceOffsetValues = [0.0]

    # init the 'Plane' selected for 'SliceType'
    Slice_1.SliceType.Origin = [X_SlicePosition, Y_SlicePosition, Z_SlicePosition]

    # init the 'Plane' selected for 'HyperTreeGridSlicer'
    Slice_1.HyperTreeGridSlicer.Origin = [X_SlicePosition, Y_SlicePosition, Z_SlicePosition]

    # show data in view
    Slice_1Display = Show(Slice_1, RenderView_1, 'GeometryRepresentation')
    
    # Hide bone slice color bar
    SliceDisplay.SetScalarBarVisibility(RenderView_1, False)    
    
    # get separate color transfer function/color map for 'MetaImage'
    MetaImageLUT = GetColorTransferFunction('MetaImage', Slice_1Display, separate=True)
    MetaImageLUT.RGBPoints = [0.18266108632087708, 0.0, 0.0, 0.5625, 0.34695134688070417, 0.0, 0.0, 1.0, 0.7224726817529351, 0.0, 1.0, 1.0, 0.9102329795355946, 0.5, 1.0, 0.5, 1.097993277318254, 1.0, 1.0, 0.0, 1.473514612190485, 1.0, 0.0, 0.0, 1.6612749099731445, 0.5, 0.0, 0.0]
    MetaImageLUT.ColorSpace = 'RGB'
    MetaImageLUT.ScalarRangeInitialized = 1.0
    MetaImageLUT.NanColor = [1.0, 0.0, 0.0]
    MetaImageLUT.NanOpacity = 0.0


    # get separate opacity transfer function/opacity map for 'MetaImage'
    MetaImagePWF = GetOpacityTransferFunction('MetaImage', Slice_1Display, separate=True)
    MetaImagePWF.Points = [0.18266108632087708, 0.0, 0.5, 0.0, 1.6612749099731445, 1.0, 0.5, 0.0]
    MetaImagePWF.ScalarRangeInitialized = 1

    # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
    MetaImageLUT.ApplyPreset('Jet', True)

    # Rescale transfer function
    MetaImageLUT.RescaleTransferFunction(Min, Max)

    # Rescale transfer function
    MetaImagePWF.RescaleTransferFunction(Min, Max)

    # get color legend/bar for MetaImageLUT in view RenderView_1
    MetaImageLUTColorBar = GetScalarBar(MetaImageLUT, RenderView_1)
    MetaImageLUTColorBar.AutoOrient = 0
    MetaImageLUTColorBar.Orientation = 'Horizontal'
    MetaImageLUTColorBar.WindowLocation = 'LowerCenter'
    MetaImageLUTColorBar.ComponentTitle = ''

    # Properties modified on MetaImageLUTColorBar
    MetaImageLUTColorBar.Title = ScalarValues_1[:-4] + ' (-)'
    MetaImageLUTColorBar.TitleFontFamily = 'Times'
    MetaImageLUTColorBar.TitleFontSize = 24
    MetaImageLUTColorBar.LabelFontFamily = 'Times'
    MetaImageLUTColorBar.LabelFontSize = 24

    # Properties modified on MetaImageLUTColorBar
    MetaImageLUTColorBar.AddRangeLabels = 0

    # trace defaults for the display properties.
    Slice_1Display.Representation = 'Surface'
    Slice_1Display.ColorArrayName = ['POINTS', 'MetaImage']
    Slice_1Display.LookupTable = MetaImageLUT
    Slice_1Display.SelectTCoordArray = 'None'
    Slice_1Display.SelectNormalArray = 'None'
    Slice_1Display.SelectTangentArray = 'None'
    Slice_1Display.OSPRayScaleArray = 'MetaImage'
    Slice_1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    Slice_1Display.SelectOrientationVectors = 'None'
    Slice_1Display.ScaleFactor = 47.400000000000006
    Slice_1Display.SelectScaleArray = 'MetaImage'
    Slice_1Display.GlyphType = 'Arrow'
    Slice_1Display.GlyphTableIndexArray = 'MetaImage'
    Slice_1Display.GaussianRadius = 2.37
    Slice_1Display.SetScaleArray = ['POINTS', 'MetaImage']
    Slice_1Display.ScaleTransferFunction = 'PiecewiseFunction'
    Slice_1Display.OpacityArray = ['POINTS', 'MetaImage']
    Slice_1Display.OpacityTransferFunction = 'PiecewiseFunction'
    Slice_1Display.DataAxesGrid = 'GridAxesRepresentation'
    Slice_1Display.PolarAxes = 'PolarAxesRepresentation'
    
    # Modify opacity
    Slice_1Display.Opacity = ImageOpacity

    # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
    Slice_1Display.OSPRayScaleFunction.Points = [-82.39119720458984, 0.0, 0.5, 0.0, 82.35810089111328, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    Slice_1Display.ScaleTransferFunction.Points = [-1699.5, 0.0, 0.5, 0.0, 14097.0, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    Slice_1Display.OpacityTransferFunction.Points = [-1699.5, 0.0, 0.5, 0.0, 14097.0, 1.0, 0.5, 0.0]

    # show color bar/color legend
    Slice_1Display.SetScalarBarVisibility(RenderView_1, True)

    # set active source
    SetActiveSource(Image_1)

    # update the view to ensure updated data information
    RenderView_1.Update()


    # get layout
    layout1 = GetLayout()

    # layout/tab size in pixels
    layout1.SetSize(998, 522)

    # save screenshot
    SaveScreenshot(SaveDir+'02_Moving_'+ScalarValues_1[:-4]+'_YZPlane.png', RenderView_1, ImageResolution=[998, 522])

    return Image_1, Slice_1, Slice_1Display



# 02 Initialize variables and environment
WorkingDirectory = '/home/mathieu/Documents/Post-Msc2/04_Results/04_Registration/'
Samples = [Folder for Folder in os.listdir(WorkingDirectory) if os.path.isdir(WorkingDirectory+Folder)]
Samples.sort()

## Get active view
RenderView_1 = GetActiveViewOrCreate('RenderView')

## Set camera placement
RenderView_1.InteractionMode = '2D'
RenderView_1.CameraPosition = [-1200, 230, 230]
RenderView_1.CameraFocalPoint = [210, 230, 200]
RenderView_1.CameraViewUp = [0, 0, -1]
RenderView_1.CameraParallelScale = 250



# 03 Print pictures of the slices
#Samples = [Samples[1],Samples[9],Samples[20]]
for Sample in Samples:
    
    print('\nPrint images for sample ' + Sample + ' ...')
    SaveDir = WorkingDirectory + Sample + '/'

    ## Fixed and moving slices
    MHDImage = 'FixedImage.mhd'
    FixedImage, FixedSlice, FixedSlice_Display, SlicePosition = PrintSlice(MHDImage,RenderView_1,SaveDir)

    Hide(FixedSlice, RenderView_1)

    MHDImage = 'RigidResult.mhd'
    MovingImage, MovingSlice, MovingSlice_Display, SlicePosition = PrintSlice(MHDImage,RenderView_1,SaveDir)



    ## Add J and F_tilde slices in transparency
    ScalarValues_1 = 'J.mhd'
    Min, Max = 0, 2
    JImage, JSlice, JSlice_Display = PrintDeformations(MovingSlice_Display,SlicePosition,ScalarValues_1,Min,Max,RenderView_1,SaveDir,ImageOpacity=0.5)

    Hide(JSlice, RenderView_1)

    ScalarValues_1 = 'F_Tilde.mhd'
    Min, Max = 1.7, 2
    FImage, FSlice, FSlice_Display = PrintDeformations(MovingSlice_Display,SlicePosition,ScalarValues_1,Min,Max,RenderView_1,SaveDir,ImageOpacity=0.5)

    ## Hide moving slice and show fixed slice
    Hide(MovingSlice, RenderView_1)
    Show(FixedSlice, RenderView_1)
    RenderView_1.Update()
    SaveScreenshot(SaveDir+'01_Fixed_F_Tilde_YZPlane.png', RenderView_1, ImageResolution=[998, 522])

    ## Hide F_Tilde and J
    Hide(FSlice, RenderView_1)
    Show(JSlice, RenderView_1)
    JSlice_Display.SetScalarBarVisibility(RenderView_1, True)
    RenderView_1.Update()
    SaveScreenshot(SaveDir+'01_Fixed_J_YZPlane.png', RenderView_1, ImageResolution=[998, 522])
    JSlice_Display.SetScalarBarVisibility(RenderView_1, False)
    RenderView_1.Update()


    ## Clean variables and go to next sample
    Delete(FImage)
    del FImage

    Delete(FSlice)
    del FSlice

    Delete(JImage)
    del JImage

    Delete(JSlice)
    del JSlice

    Delete(MovingImage)
    del MovingImage

    Delete(MovingSlice)
    del MovingSlice

    Delete(FixedImage)
    del FixedImage

    Delete(FixedSlice)
    del FixedSlice
