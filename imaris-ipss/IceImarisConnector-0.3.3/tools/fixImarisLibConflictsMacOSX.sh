#!/bin/bash

# When starting Imaris 7.6 from MATLAB on Mac OS X, Imaris crashes
# because it tries to load the Qt library from MATLAB which are at a
# different version. Running this script once, changes the Qt library
# names and pointers to workaround the conflict.
#
# Aaron Ponti, 2012/11/13

ImarisPath="/Applications/Imaris 7.6.0.app"

cd "/"

# Rename libraries
sudo mv "${ImarisPath}/Contents/Frameworks/QtCore.framework/Versions/4/QtCore" \
	"${ImarisPath}/Contents/Frameworks/QtCore.framework/Versions/4/QtCore48"
sudo mv "${ImarisPath}/Contents/Frameworks/QtDeclarative.framework/Versions/4/QtDeclarative" \
	"${ImarisPath}/Contents/Frameworks/QtDeclarative.framework/Versions/4/QtDeclarative48"
sudo mv "${ImarisPath}/Contents/Frameworks/QtGui.framework/Versions/4/QtGui" \
	"${ImarisPath}/Contents/Frameworks/QtGui.framework/Versions/4/QtGui48"
sudo mv "${ImarisPath}/Contents/Frameworks/QtNetwork.framework/Versions/4/QtNetwork" \
	"${ImarisPath}/Contents/Frameworks/QtNetwork.framework/Versions/4/QtNetwork48"
sudo mv "${ImarisPath}/Contents/Frameworks/QtOpenGL.framework/Versions/4/QtOpenGL" \
	"${ImarisPath}/Contents/Frameworks/QtOpenGL.framework/Versions/4/QtOpenGL48"
sudo mv "${ImarisPath}/Contents/Frameworks/QtScript.framework/Versions/4/QtScript" \
	"${ImarisPath}/Contents/Frameworks/QtScript.framework/Versions/4/QtScript48"
sudo mv "${ImarisPath}/Contents/Frameworks/QtSql.framework/Versions/4/QtSql" \
	"${ImarisPath}/Contents/Frameworks/QtSql.framework/Versions/4/QtSql48"
sudo mv "${ImarisPath}/Contents/Frameworks/QtSvg.framework/Versions/4/QtSvg" \
	"${ImarisPath}/Contents/Frameworks/QtSvg.framework/Versions/4/QtSvg48"
sudo mv "${ImarisPath}/Contents/Frameworks/QtXml.framework/Versions/4/QtXml" \
	"${ImarisPath}/Contents/Frameworks/QtXml.framework/Versions/4/QtXml48"
sudo mv "${ImarisPath}/Contents/Frameworks/QtXmlPatterns.framework/Versions/4/QtXmlPatterns" \
	"${ImarisPath}/Contents/Frameworks/QtXmlPatterns.framework/Versions/4/QtXmlPatterns48"

# Change ID
sudo install_name_tool -id @loader_path/QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/Frameworks/QtCore.framework/Versions/4/QtCore48"	
sudo install_name_tool -id @loader_path/QtDeclarative.framework/Versions/4/QtDeclarative48 \
	"${ImarisPath}/Contents/Frameworks/QtDeclarative.framework/Versions/4/QtDeclarative48"
sudo install_name_tool -id @loader_path/QtGui.framework/Versions/4/QtGui48 \
	"${ImarisPath}/Contents/Frameworks/QtGui.framework/Versions/4/QtGui48"
sudo install_name_tool -id @loader_path/QtNetwork.framework/Versions/4/QtNetwork48 \
	"${ImarisPath}/Contents/Frameworks/QtNetwork.framework/Versions/4/QtNetwork48"
sudo install_name_tool -id @loader_path/QtOpenGL.framework/Versions/4/QtOpenGL48 \
	"${ImarisPath}/Contents/Frameworks/QtOpenGL.framework/Versions/4/QtOpenGL48"
sudo install_name_tool -id @loader_path/QtScript.framework/Versions/4/QtScript48 \
	"${ImarisPath}/Contents/Frameworks/QtScript.framework/Versions/4/QtScript48"
sudo install_name_tool -id @loader_path/QtSql.framework/Versions/4/QtSql48 \
	"${ImarisPath}/Contents/Frameworks/QtSql.framework/Versions/4/QtSql48"
sudo install_name_tool -id @loader_path/QtSvg.framework/Versions/4/QtSvg48 \
	"${ImarisPath}/Contents/Frameworks/QtSvg.framework/Versions/4/QtSvg48"
sudo install_name_tool -id @loader_path/QtXml.framework/Versions/4/QtXml48 \
	"${ImarisPath}/Contents/Frameworks/QtXml.framework/Versions/4/QtXml48"
sudo install_name_tool -id @loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns48 \
	"${ImarisPath}/Contents/Frameworks/QtXmlPatterns.framework/Versions/4/QtXmlPatterns48"

# Change the references of the Qt frameworks with each other
sudo install_name_tool -change \
	@loader_path/../../../QtCore.framework/Versions/4/QtCore \
	@loader_path/../../../QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/MacOS/../Frameworks/QtGui.framework/Versions/4/QtGui48"
sudo install_name_tool -change \
	@loader_path/../../../QtCore.framework/Versions/4/QtCore \
	@loader_path/../../../QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/MacOS/../Frameworks/QtNetwork.framework/Versions/4/QtNetwork48"
sudo install_name_tool -change \
	@loader_path/../../../QtGui.framework/Versions/4/QtGui \
	@loader_path/../../../QtGui.framework/Versions/4/QtGui48 \
	"${ImarisPath}/Contents/MacOS/../Frameworks/QtOpenGL.framework/Versions/4/QtOpenGL48"
sudo install_name_tool -change \
	@loader_path/../../../QtCore.framework/Versions/4/QtCore \
	@loader_path/../../../QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/MacOS/../Frameworks/QtOpenGL.framework/Versions/4/QtOpenGL48"
sudo install_name_tool -change \
	@loader_path/../../../QtCore.framework/Versions/4/QtCore \
	@loader_path/../../../QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/MacOS/../Frameworks/QtXML.framework/Versions/4/QtXML48"
sudo install_name_tool -change \
	@loader_path/../../../QtCore.framework/Versions/4/QtCore \
	@loader_path/../../../QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/MacOS/../Frameworks/QtSql.framework/Versions/4/QtSql48"

# Change the references of other libs
sudo install_name_tool -change @loader_path/QtCore.framework/Versions/4/QtCore \
	@loader_path/QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/Frameworks/libCg.2.1.dylib"
sudo install_name_tool -change @loader_path/QtCore.framework/Versions/4/QtCore \
	@loader_path/QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/Frameworks/libCoin.60.dylib"
sudo install_name_tool -change @loader_path/QtCore.framework/Versions/4/QtCore \
	@loader_path/QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/Frameworks/libIce.34.dylib"
sudo install_name_tool -change @loader_path/QtCore.framework/Versions/4/QtCore \
	@loader_path/QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/Frameworks/libSBReadFile.dylib"
sudo install_name_tool -change @loader_path/QtCore.framework/Versions/4/QtCore \
	@loader_path/QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/Frameworks/libSoQt.20.dylib"
sudo install_name_tool -change @loader_path/QtCore.framework/Versions/4/QtCore \
	@loader_path/QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/Frameworks/libavcodec.dylib"
sudo install_name_tool -change @loader_path/QtCore.framework/Versions/4/QtCore \
	@loader_path/QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/Frameworks/libavformat.dylib"
sudo install_name_tool -change @loader_path/QtCore.framework/Versions/4/QtCore \
	@loader_path/QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/Frameworks/libavutil.dylib"
sudo install_name_tool -change @loader_path/QtCore.framework/Versions/4/QtCore \
	@loader_path/QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/Frameworks/libbpBatchIcePlugin.7.6.dylib"
sudo install_name_tool -change @loader_path/QtCore.framework/Versions/4/QtCore \
	@loader_path/QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/Frameworks/libbpImarisBatchAgent.7.6.dylib"
sudo install_name_tool -change @loader_path/QtCore.framework/Versions/4/QtCore \
	@loader_path/QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/Frameworks/libbpcore.7.6.dylib"
sudo install_name_tool -change @loader_path/QtCore.framework/Versions/4/QtCore \
	@loader_path/QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/Frameworks/libbpexception.7.6.dylib"
sudo install_name_tool -change @loader_path/QtCore.framework/Versions/4/QtCore \
	@loader_path/QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/Frameworks/libbpfileio.7.6.dylib"
sudo install_name_tool -change @loader_path/QtCore.framework/Versions/4/QtCore \
	@loader_path/QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/Frameworks/libexpat.1.dylib"
sudo install_name_tool -change @loader_path/QtCore.framework/Versions/4/QtCore \
	@loader_path/QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/Frameworks/libfreeimage-3.15.0.dylib"
sudo install_name_tool -change @loader_path/QtCore.framework/Versions/4/QtCore \
	@loader_path/QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/Frameworks/libhdf5.7.dylib"
sudo install_name_tool -change @loader_path/QtCore.framework/Versions/4/QtCore \
	@loader_path/QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/Frameworks/libswscale.dylib"
###
sudo install_name_tool -change @loader_path/QtDeclarative.framework/Versions/4/QtDeclarative \
	@loader_path/QtDeclarative.framework/Versions/4/QtDeclarative48 \
	"${ImarisPath}/Contents/Frameworks/libCg.2.1.dylib"
sudo install_name_tool -change @loader_path/QtDeclarative.framework/Versions/4/QtDeclarative \
	@loader_path/QtDeclarative.framework/Versions/4/QtDeclarative48 \
	"${ImarisPath}/Contents/Frameworks/libCoin.60.dylib"
sudo install_name_tool -change @loader_path/QtDeclarative.framework/Versions/4/QtDeclarative \
	@loader_path/QtDeclarative.framework/Versions/4/QtDeclarative48 \
	"${ImarisPath}/Contents/Frameworks/libIce.34.dylib"
sudo install_name_tool -change @loader_path/QtDeclarative.framework/Versions/4/QtDeclarative \
	@loader_path/QtDeclarative.framework/Versions/4/QtDeclarative48 \
	"${ImarisPath}/Contents/Frameworks/libSBReadFile.dylib"
sudo install_name_tool -change @loader_path/QtDeclarative.framework/Versions/4/QtDeclarative \
	@loader_path/QtDeclarative.framework/Versions/4/QtDeclarative48 \
	"${ImarisPath}/Contents/Frameworks/libSoQt.20.dylib"
sudo install_name_tool -change @loader_path/QtDeclarative.framework/Versions/4/QtDeclarative \
	@loader_path/QtDeclarative.framework/Versions/4/QtDeclarative48 \
	"${ImarisPath}/Contents/Frameworks/libavcodec.dylib"
sudo install_name_tool -change @loader_path/QtDeclarative.framework/Versions/4/QtDeclarative \
	@loader_path/QtDeclarative.framework/Versions/4/QtDeclarative48 \
	"${ImarisPath}/Contents/Frameworks/libavformat.dylib"
sudo install_name_tool -change @loader_path/QtDeclarative.framework/Versions/4/QtDeclarative \
	@loader_path/QtDeclarative.framework/Versions/4/QtDeclarative48 \
	"${ImarisPath}/Contents/Frameworks/libavutil.dylib"
sudo install_name_tool -change @loader_path/QtDeclarative.framework/Versions/4/QtDeclarative \
	@loader_path/QtDeclarative.framework/Versions/4/QtDeclarative48 \
	"${ImarisPath}/Contents/Frameworks/libbpBatchIcePlugin.7.6.dylib"
sudo install_name_tool -change @loader_path/QtDeclarative.framework/Versions/4/QtDeclarative \
	@loader_path/QtDeclarative.framework/Versions/4/QtDeclarative48 \
	"${ImarisPath}/Contents/Frameworks/libbpImarisBatchAgent.7.6.dylib"
sudo install_name_tool -change @loader_path/QtDeclarative.framework/Versions/4/QtDeclarative \
	@loader_path/QtDeclarative.framework/Versions/4/QtDeclarative48 \
	"${ImarisPath}/Contents/Frameworks/libbpcore.7.6.dylib"
sudo install_name_tool -change @loader_path/QtDeclarative.framework/Versions/4/QtDeclarative \
	@loader_path/QtDeclarative.framework/Versions/4/QtDeclarative48 \
	"${ImarisPath}/Contents/Frameworks/libbpexception.7.6.dylib"
sudo install_name_tool -change @loader_path/QtDeclarative.framework/Versions/4/QtDeclarative \
	@loader_path/QtDeclarative.framework/Versions/4/QtDeclarative48 \
	"${ImarisPath}/Contents/Frameworks/libbpfileio.7.6.dylib"
sudo install_name_tool -change @loader_path/QtDeclarative.framework/Versions/4/QtDeclarative \
	@loader_path/QtDeclarative.framework/Versions/4/QtDeclarative48 \
	"${ImarisPath}/Contents/Frameworks/libexpat.1.dylib"
sudo install_name_tool -change @loader_path/QtDeclarative.framework/Versions/4/QtDeclarative \
	@loader_path/QtDeclarative.framework/Versions/4/QtDeclarative48 \
	"${ImarisPath}/Contents/Frameworks/libfreeimage-3.15.0.dylib"
sudo install_name_tool -change @loader_path/QtDeclarative.framework/Versions/4/QtDeclarative \
	@loader_path/QtDeclarative.framework/Versions/4/QtDeclarative48 \
	"${ImarisPath}/Contents/Frameworks/libhdf5.7.dylib"
sudo install_name_tool -change @loader_path/QtDeclarative.framework/Versions/4/QtDeclarative \
	@loader_path/QtDeclarative.framework/Versions/4/QtDeclarative48 \
	"${ImarisPath}/Contents/Frameworks/libswscale.dylib"
###
sudo install_name_tool -change @loader_path/QtGui.framework/Versions/4/QtGui \
	@loader_path/QtGui.framework/Versions/4/QtGui48 \
	"${ImarisPath}/Contents/Frameworks/libCg.2.1.dylib"
sudo install_name_tool -change @loader_path/QtGui.framework/Versions/4/QtGui \
	@loader_path/QtGui.framework/Versions/4/QtGui48 \
	"${ImarisPath}/Contents/Frameworks/libCoin.60.dylib"
sudo install_name_tool -change @loader_path/QtGui.framework/Versions/4/QtGui \
	@loader_path/QtGui.framework/Versions/4/QtGui48 \
	"${ImarisPath}/Contents/Frameworks/libIce.34.dylib"
sudo install_name_tool -change @loader_path/QtGui.framework/Versions/4/QtGui \
	@loader_path/QtGui.framework/Versions/4/QtGui48 \
	"${ImarisPath}/Contents/Frameworks/libSBReadFile.dylib"
sudo install_name_tool -change @loader_path/QtGui.framework/Versions/4/QtGui \
	@loader_path/QtGui.framework/Versions/4/QtGui48 \
	"${ImarisPath}/Contents/Frameworks/libSoQt.20.dylib"
sudo install_name_tool -change @loader_path/QtGui.framework/Versions/4/QtGui \
	@loader_path/QtGui.framework/Versions/4/QtGui48 \
	"${ImarisPath}/Contents/Frameworks/libavcodec.dylib"
sudo install_name_tool -change @loader_path/QtGui.framework/Versions/4/QtGui \
	@loader_path/QtGui.framework/Versions/4/QtGui48 \
	"${ImarisPath}/Contents/Frameworks/libavformat.dylib"
sudo install_name_tool -change @loader_path/QtGui.framework/Versions/4/QtGui \
	@loader_path/QtGui.framework/Versions/4/QtGui48 \
	"${ImarisPath}/Contents/Frameworks/libavutil.dylib"
sudo install_name_tool -change @loader_path/QtGui.framework/Versions/4/QtGui \
	@loader_path/QtGui.framework/Versions/4/QtGui48 \
	"${ImarisPath}/Contents/Frameworks/libbpBatchIcePlugin.7.6.dylib"
sudo install_name_tool -change @loader_path/QtGui.framework/Versions/4/QtGui \
	@loader_path/QtGui.framework/Versions/4/QtGui48 \
	"${ImarisPath}/Contents/Frameworks/libbpImarisBatchAgent.7.6.dylib"
sudo install_name_tool -change @loader_path/QtGui.framework/Versions/4/QtGui \
	@loader_path/QtGui.framework/Versions/4/QtGui48 \
	"${ImarisPath}/Contents/Frameworks/libbpcore.7.6.dylib"
sudo install_name_tool -change @loader_path/QtGui.framework/Versions/4/QtGui \
	@loader_path/QtGui.framework/Versions/4/QtGui48 \
	"${ImarisPath}/Contents/Frameworks/libbpexception.7.6.dylib"
sudo install_name_tool -change @loader_path/QtGui.framework/Versions/4/QtGui \
	@loader_path/QtGui.framework/Versions/4/QtGui48 \
	"${ImarisPath}/Contents/Frameworks/libbpfileio.7.6.dylib"
sudo install_name_tool -change @loader_path/QtGui.framework/Versions/4/QtGui \
	@loader_path/QtGui.framework/Versions/4/QtGui48 \
	"${ImarisPath}/Contents/Frameworks/libexpat.1.dylib"
sudo install_name_tool -change @loader_path/QtGui.framework/Versions/4/QtGui \
	@loader_path/QtGui.framework/Versions/4/QtGui48 \
	"${ImarisPath}/Contents/Frameworks/libfreeimage-3.15.0.dylib"
sudo install_name_tool -change @loader_path/QtGui.framework/Versions/4/QtGui \
	@loader_path/QtGui.framework/Versions/4/QtGui48 \
	"${ImarisPath}/Contents/Frameworks/libhdf5.7.dylib"
sudo install_name_tool -change @loader_path/QtGui.framework/Versions/4/QtGui \
	@loader_path/QtGui.framework/Versions/4/QtGui48 \
	"${ImarisPath}/Contents/Frameworks/libswscale.dylib"
###
sudo install_name_tool -change @loader_path/QtNetwork.framework/Versions/4/QtNetwork \
	@loader_path/QtNetwork.framework/Versions/4/QtNetwork48 \
	"${ImarisPath}/Contents/Frameworks/libCg.2.1.dylib"
sudo install_name_tool -change @loader_path/QtNetwork.framework/Versions/4/QtNetwork \
	@loader_path/QtNetwork.framework/Versions/4/QtNetwork48 \
	"${ImarisPath}/Contents/Frameworks/libCoin.60.dylib"
sudo install_name_tool -change @loader_path/QtNetwork.framework/Versions/4/QtNetwork \
	@loader_path/QtNetwork.framework/Versions/4/QtNetwork48 \
	"${ImarisPath}/Contents/Frameworks/libIce.34.dylib"
sudo install_name_tool -change @loader_path/QtNetwork.framework/Versions/4/QtNetwork \
	@loader_path/QtNetwork.framework/Versions/4/QtNetwork48 \
	"${ImarisPath}/Contents/Frameworks/libSBReadFile.dylib"
sudo install_name_tool -change @loader_path/QtNetwork.framework/Versions/4/QtNetwork \
	@loader_path/QtNetwork.framework/Versions/4/QtNetwork48 \
	"${ImarisPath}/Contents/Frameworks/libSoQt.20.dylib"
sudo install_name_tool -change @loader_path/QtNetwork.framework/Versions/4/QtNetwork \
	@loader_path/QtNetwork.framework/Versions/4/QtNetwork48 \
	"${ImarisPath}/Contents/Frameworks/libavcodec.dylib"
sudo install_name_tool -change @loader_path/QtNetwork.framework/Versions/4/QtNetwork \
	@loader_path/QtNetwork.framework/Versions/4/QtNetwork48 \
	"${ImarisPath}/Contents/Frameworks/libavformat.dylib"
sudo install_name_tool -change @loader_path/QtNetwork.framework/Versions/4/QtNetwork \
	@loader_path/QtNetwork.framework/Versions/4/QtNetwork48 \
	"${ImarisPath}/Contents/Frameworks/libavutil.dylib"
sudo install_name_tool -change @loader_path/QtNetwork.framework/Versions/4/QtNetwork \
	@loader_path/QtNetwork.framework/Versions/4/QtNetwork48 \
	"${ImarisPath}/Contents/Frameworks/libbpBatchIcePlugin.7.6.dylib"
sudo install_name_tool -change @loader_path/QtNetwork.framework/Versions/4/QtNetwork \
	@loader_path/QtNetwork.framework/Versions/4/QtNetwork48 \
	"${ImarisPath}/Contents/Frameworks/libbpImarisBatchAgent.7.6.dylib"
sudo install_name_tool -change @loader_path/QtNetwork.framework/Versions/4/QtNetwork \
	@loader_path/QtNetwork.framework/Versions/4/QtNetwork48 \
	"${ImarisPath}/Contents/Frameworks/libbpcore.7.6.dylib"
sudo install_name_tool -change @loader_path/QtNetwork.framework/Versions/4/QtNetwork \
	@loader_path/QtNetwork.framework/Versions/4/QtNetwork48 \
	"${ImarisPath}/Contents/Frameworks/libbpexception.7.6.dylib"
sudo install_name_tool -change @loader_path/QtNetwork.framework/Versions/4/QtNetwork \
	@loader_path/QtNetwork.framework/Versions/4/QtNetwork48 \
	"${ImarisPath}/Contents/Frameworks/libbpfileio.7.6.dylib"
sudo install_name_tool -change @loader_path/QtNetwork.framework/Versions/4/QtNetwork \
	@loader_path/QtNetwork.framework/Versions/4/QtNetwork48 \
	"${ImarisPath}/Contents/Frameworks/libexpat.1.dylib"
sudo install_name_tool -change @loader_path/QtNetwork.framework/Versions/4/QtNetwork \
	@loader_path/QtNetwork.framework/Versions/4/QtNetwork48 \
	"${ImarisPath}/Contents/Frameworks/libfreeimage-3.15.0.dylib"
sudo install_name_tool -change @loader_path/QtNetwork.framework/Versions/4/QtNetwork \
	@loader_path/QtNetwork.framework/Versions/4/QtNetwork48 \
	"${ImarisPath}/Contents/Frameworks/libhdf5.7.dylib"
sudo install_name_tool -change @loader_path/QtNetwork.framework/Versions/4/QtNetwork \
	@loader_path/QtNetwork.framework/Versions/4/QtNetwork48 \
	"${ImarisPath}/Contents/Frameworks/libswscale.dylib"
###
sudo install_name_tool -change @loader_path/QtOpenGL.framework/Versions/4/QtOpenGL \
	@loader_path/QtOpenGL.framework/Versions/4/QtOpenGL48 \
	"${ImarisPath}/Contents/Frameworks/libCg.2.1.dylib"
sudo install_name_tool -change @loader_path/QtOpenGL.framework/Versions/4/QtOpenGL \
	@loader_path/QtOpenGL.framework/Versions/4/QtOpenGL48 \
	"${ImarisPath}/Contents/Frameworks/libCoin.60.dylib"
sudo install_name_tool -change @loader_path/QtOpenGL.framework/Versions/4/QtOpenGL \
	@loader_path/QtOpenGL.framework/Versions/4/QtOpenGL48 \
	"${ImarisPath}/Contents/Frameworks/libIce.34.dylib"
sudo install_name_tool -change @loader_path/QtOpenGL.framework/Versions/4/QtOpenGL \
	@loader_path/QtOpenGL.framework/Versions/4/QtOpenGL48 \
	"${ImarisPath}/Contents/Frameworks/libSBReadFile.dylib"
sudo install_name_tool -change @loader_path/QtOpenGL.framework/Versions/4/QtOpenGL \
	@loader_path/QtOpenGL.framework/Versions/4/QtOpenGL48 \
	"${ImarisPath}/Contents/Frameworks/libSoQt.20.dylib"
sudo install_name_tool -change @loader_path/QtOpenGL.framework/Versions/4/QtOpenGL \
	@loader_path/QtOpenGL.framework/Versions/4/QtOpenGL48 \
	"${ImarisPath}/Contents/Frameworks/libavcodec.dylib"
sudo install_name_tool -change @loader_path/QtOpenGL.framework/Versions/4/QtOpenGL \
	@loader_path/QtOpenGL.framework/Versions/4/QtOpenGL48 \
	"${ImarisPath}/Contents/Frameworks/libavformat.dylib"
sudo install_name_tool -change @loader_path/QtOpenGL.framework/Versions/4/QtOpenGL \
	@loader_path/QtOpenGL.framework/Versions/4/QtOpenGL48 \
	"${ImarisPath}/Contents/Frameworks/libavutil.dylib"
sudo install_name_tool -change @loader_path/QtOpenGL.framework/Versions/4/QtOpenGL \
	@loader_path/QtOpenGL.framework/Versions/4/QtOpenGL48 \
	"${ImarisPath}/Contents/Frameworks/libbpBatchIcePlugin.7.6.dylib"
sudo install_name_tool -change @loader_path/QtOpenGL.framework/Versions/4/QtOpenGL \
	@loader_path/QtOpenGL.framework/Versions/4/QtOpenGL48 \
	"${ImarisPath}/Contents/Frameworks/libbpImarisBatchAgent.7.6.dylib"
sudo install_name_tool -change @loader_path/QtOpenGL.framework/Versions/4/QtOpenGL \
	@loader_path/QtOpenGL.framework/Versions/4/QtOpenGL48 \
	"${ImarisPath}/Contents/Frameworks/libbpcore.7.6.dylib"
sudo install_name_tool -change @loader_path/QtOpenGL.framework/Versions/4/QtOpenGL \
	@loader_path/QtOpenGL.framework/Versions/4/QtOpenGL48 \
	"${ImarisPath}/Contents/Frameworks/libbpexception.7.6.dylib"
sudo install_name_tool -change @loader_path/QtOpenGL.framework/Versions/4/QtOpenGL \
	@loader_path/QtOpenGL.framework/Versions/4/QtOpenGL48 \
	"${ImarisPath}/Contents/Frameworks/libbpfileio.7.6.dylib"
sudo install_name_tool -change @loader_path/QtOpenGL.framework/Versions/4/QtOpenGL \
	@loader_path/QtOpenGL.framework/Versions/4/QtOpenGL48 \
	"${ImarisPath}/Contents/Frameworks/libexpat.1.dylib"
sudo install_name_tool -change @loader_path/QtOpenGL.framework/Versions/4/QtOpenGL \
	@loader_path/QtOpenGL.framework/Versions/4/QtOpenGL48 \
	"${ImarisPath}/Contents/Frameworks/libfreeimage-3.15.0.dylib"
sudo install_name_tool -change @loader_path/QtOpenGL.framework/Versions/4/QtOpenGL \
	@loader_path/QtOpenGL.framework/Versions/4/QtOpenGL48 \
	"${ImarisPath}/Contents/Frameworks/libhdf5.7.dylib"
sudo install_name_tool -change @loader_path/QtOpenGL.framework/Versions/4/QtOpenGL \
	@loader_path/QtOpenGL.framework/Versions/4/QtOpenGL48 \
	"${ImarisPath}/Contents/Frameworks/libswscale.dylib"
###
sudo install_name_tool -change @loader_path/QtScript.framework/Versions/4/QtScript \
	@loader_path/QtScript.framework/Versions/4/QtScript48 \
	"${ImarisPath}/Contents/Frameworks/libCg.2.1.dylib"
sudo install_name_tool -change @loader_path/QtScript.framework/Versions/4/QtScript \
	@loader_path/QtScript.framework/Versions/4/QtScript48 \
	"${ImarisPath}/Contents/Frameworks/libCoin.60.dylib"
sudo install_name_tool -change @loader_path/QtScript.framework/Versions/4/QtScript \
	@loader_path/QtScript.framework/Versions/4/QtScript48 \
	"${ImarisPath}/Contents/Frameworks/libIce.34.dylib"
sudo install_name_tool -change @loader_path/QtScript.framework/Versions/4/QtScript \
	@loader_path/QtScript.framework/Versions/4/QtScript48 \
	"${ImarisPath}/Contents/Frameworks/libSBReadFile.dylib"
sudo install_name_tool -change @loader_path/QtScript.framework/Versions/4/QtScript \
	@loader_path/QtScript.framework/Versions/4/QtScript48 \
	"${ImarisPath}/Contents/Frameworks/libSoQt.20.dylib"
sudo install_name_tool -change @loader_path/QtScript.framework/Versions/4/QtScript \
	@loader_path/QtScript.framework/Versions/4/QtScript48 \
	"${ImarisPath}/Contents/Frameworks/libavcodec.dylib"
sudo install_name_tool -change @loader_path/QtScript.framework/Versions/4/QtScript \
	@loader_path/QtScript.framework/Versions/4/QtScript48 \
	"${ImarisPath}/Contents/Frameworks/libavformat.dylib"
sudo install_name_tool -change @loader_path/QtScript.framework/Versions/4/QtScript \
	@loader_path/QtScript.framework/Versions/4/QtScript48 \
	"${ImarisPath}/Contents/Frameworks/libavutil.dylib"
sudo install_name_tool -change @loader_path/QtScript.framework/Versions/4/QtScript \
	@loader_path/QtScript.framework/Versions/4/QtScript48 \
	"${ImarisPath}/Contents/Frameworks/libbpBatchIcePlugin.7.6.dylib"
sudo install_name_tool -change @loader_path/QtScript.framework/Versions/4/QtScript \
	@loader_path/QtScript.framework/Versions/4/QtScript48 \
	"${ImarisPath}/Contents/Frameworks/libbpImarisBatchAgent.7.6.dylib"
sudo install_name_tool -change @loader_path/QtScript.framework/Versions/4/QtScript \
	@loader_path/QtScript.framework/Versions/4/QtScript48 \
	"${ImarisPath}/Contents/Frameworks/libbpcore.7.6.dylib"
sudo install_name_tool -change @loader_path/QtScript.framework/Versions/4/QtScript \
	@loader_path/QtScript.framework/Versions/4/QtScript48 \
	"${ImarisPath}/Contents/Frameworks/libbpexception.7.6.dylib"
sudo install_name_tool -change @loader_path/QtScript.framework/Versions/4/QtScript \
	@loader_path/QtScript.framework/Versions/4/QtScript48 \
	"${ImarisPath}/Contents/Frameworks/libbpfileio.7.6.dylib"
sudo install_name_tool -change @loader_path/QtScript.framework/Versions/4/QtScript \
	@loader_path/QtScript.framework/Versions/4/QtScript48 \
	"${ImarisPath}/Contents/Frameworks/libexpat.1.dylib"
sudo install_name_tool -change @loader_path/QtScript.framework/Versions/4/QtScript \
	@loader_path/QtScript.framework/Versions/4/QtScript48 \
	"${ImarisPath}/Contents/Frameworks/libfreeimage-3.15.0.dylib"
sudo install_name_tool -change @loader_path/QtScript.framework/Versions/4/QtScript \
	@loader_path/QtScript.framework/Versions/4/QtScript48 \
	"${ImarisPath}/Contents/Frameworks/libhdf5.7.dylib"
sudo install_name_tool -change @loader_path/QtScript.framework/Versions/4/QtScript \
	@loader_path/QtScript.framework/Versions/4/QtScript48 \
	"${ImarisPath}/Contents/Frameworks/libswscale.dylib"
###
sudo install_name_tool -change @loader_path/QtSql.framework/Versions/4/QtSql \
	@loader_path/QtSql.framework/Versions/4/QtSql48 \
	"${ImarisPath}/Contents/Frameworks/libCg.2.1.dylib"
sudo install_name_tool -change @loader_path/QtSql.framework/Versions/4/QtSql \
	@loader_path/QtSql.framework/Versions/4/QtSql48 \
	"${ImarisPath}/Contents/Frameworks/libCoin.60.dylib"
sudo install_name_tool -change @loader_path/QtSql.framework/Versions/4/QtSql \
	@loader_path/QtSql.framework/Versions/4/QtSql48 \
	"${ImarisPath}/Contents/Frameworks/libIce.34.dylib"
sudo install_name_tool -change @loader_path/QtSql.framework/Versions/4/QtSql \
	@loader_path/QtSql.framework/Versions/4/QtSql48 \
	"${ImarisPath}/Contents/Frameworks/libSBReadFile.dylib"
sudo install_name_tool -change @loader_path/QtSql.framework/Versions/4/QtSql \
	@loader_path/QtSql.framework/Versions/4/QtSql48 \
	"${ImarisPath}/Contents/Frameworks/libSoQt.20.dylib"
sudo install_name_tool -change @loader_path/QtSql.framework/Versions/4/QtSql \
	@loader_path/QtSql.framework/Versions/4/QtSql48 \
	"${ImarisPath}/Contents/Frameworks/libavcodec.dylib"
sudo install_name_tool -change @loader_path/QtSql.framework/Versions/4/QtSql \
	@loader_path/QtSql.framework/Versions/4/QtSql48 \
	"${ImarisPath}/Contents/Frameworks/libavformat.dylib"
sudo install_name_tool -change @loader_path/QtSql.framework/Versions/4/QtSql \
	@loader_path/QtSql.framework/Versions/4/QtSql48 \
	"${ImarisPath}/Contents/Frameworks/libavutil.dylib"
sudo install_name_tool -change @loader_path/QtSql.framework/Versions/4/QtSql \
	@loader_path/QtSql.framework/Versions/4/QtSql48 \
	"${ImarisPath}/Contents/Frameworks/libbpBatchIcePlugin.7.6.dylib"
sudo install_name_tool -change @loader_path/QtSql.framework/Versions/4/QtSql \
	@loader_path/QtSql.framework/Versions/4/QtSql48 \
	"${ImarisPath}/Contents/Frameworks/libbpImarisBatchAgent.7.6.dylib"
sudo install_name_tool -change @loader_path/QtSql.framework/Versions/4/QtSql \
	@loader_path/QtSql.framework/Versions/4/QtSql48 \
	"${ImarisPath}/Contents/Frameworks/libbpcore.7.6.dylib"
sudo install_name_tool -change @loader_path/QtSql.framework/Versions/4/QtSql \
	@loader_path/QtSql.framework/Versions/4/QtSql48 \
	"${ImarisPath}/Contents/Frameworks/libbpexception.7.6.dylib"
sudo install_name_tool -change @loader_path/QtSql.framework/Versions/4/QtSql \
	@loader_path/QtSql.framework/Versions/4/QtSql48 \
	"${ImarisPath}/Contents/Frameworks/libbpfileio.7.6.dylib"
sudo install_name_tool -change @loader_path/QtSql.framework/Versions/4/QtSql \
	@loader_path/QtSql.framework/Versions/4/QtSql48 \
	"${ImarisPath}/Contents/Frameworks/libexpat.1.dylib"
sudo install_name_tool -change @loader_path/QtSql.framework/Versions/4/QtSql \
	@loader_path/QtSql.framework/Versions/4/QtSql48 \
	"${ImarisPath}/Contents/Frameworks/libfreeimage-3.15.0.dylib"
sudo install_name_tool -change @loader_path/QtSql.framework/Versions/4/QtSql \
	@loader_path/QtSql.framework/Versions/4/QtSql48 \
	"${ImarisPath}/Contents/Frameworks/libhdf5.7.dylib"
sudo install_name_tool -change @loader_path/QtSql.framework/Versions/4/QtSql \
	@loader_path/QtSql.framework/Versions/4/QtSql48 \
	"${ImarisPath}/Contents/Frameworks/libswscale.dylib"
###
sudo install_name_tool -change @loader_path/QtSvg.framework/Versions/4/QtSvg \
	@loader_path/QtSvg.framework/Versions/4/QtSvg48 \
	"${ImarisPath}/Contents/Frameworks/libCg.2.1.dylib"
sudo install_name_tool -change @loader_path/QtSvg.framework/Versions/4/QtSvg \
	@loader_path/QtSvg.framework/Versions/4/QtSvg48 \
	"${ImarisPath}/Contents/Frameworks/libCoin.60.dylib"
sudo install_name_tool -change @loader_path/QtSvg.framework/Versions/4/QtSvg \
	@loader_path/QtSvg.framework/Versions/4/QtSvg48 \
	"${ImarisPath}/Contents/Frameworks/libIce.34.dylib"
sudo install_name_tool -change @loader_path/QtSvg.framework/Versions/4/QtSvg \
	@loader_path/QtSvg.framework/Versions/4/QtSvg48 \
	"${ImarisPath}/Contents/Frameworks/libSBReadFile.dylib"
sudo install_name_tool -change @loader_path/QtSvg.framework/Versions/4/QtSvg \
	@loader_path/QtSvg.framework/Versions/4/QtSvg48 \
	"${ImarisPath}/Contents/Frameworks/libSoQt.20.dylib"
sudo install_name_tool -change @loader_path/QtSvg.framework/Versions/4/QtSvg \
	@loader_path/QtSvg.framework/Versions/4/QtSvg48 \
	"${ImarisPath}/Contents/Frameworks/libavcodec.dylib"
sudo install_name_tool -change @loader_path/QtSvg.framework/Versions/4/QtSvg \
	@loader_path/QtSvg.framework/Versions/4/QtSvg48 \
	"${ImarisPath}/Contents/Frameworks/libavformat.dylib"
sudo install_name_tool -change @loader_path/QtSvg.framework/Versions/4/QtSvg \
	@loader_path/QtSvg.framework/Versions/4/QtSvg48 \
	"${ImarisPath}/Contents/Frameworks/libavutil.dylib"
sudo install_name_tool -change @loader_path/QtSvg.framework/Versions/4/QtSvg \
	@loader_path/QtSvg.framework/Versions/4/QtSvg48 \
	"${ImarisPath}/Contents/Frameworks/libbpBatchIcePlugin.7.6.dylib"
sudo install_name_tool -change @loader_path/QtSvg.framework/Versions/4/QtSvg \
	@loader_path/QtSvg.framework/Versions/4/QtSvg48 \
	"${ImarisPath}/Contents/Frameworks/libbpImarisBatchAgent.7.6.dylib"
sudo install_name_tool -change @loader_path/QtSvg.framework/Versions/4/QtSvg \
	@loader_path/QtSvg.framework/Versions/4/QtSvg48 \
	"${ImarisPath}/Contents/Frameworks/libbpcore.7.6.dylib"
sudo install_name_tool -change @loader_path/QtSvg.framework/Versions/4/QtSvg \
	@loader_path/QtSvg.framework/Versions/4/QtSvg48 \
	"${ImarisPath}/Contents/Frameworks/libbpexception.7.6.dylib"
sudo install_name_tool -change @loader_path/QtSvg.framework/Versions/4/QtSvg \
	@loader_path/QtSvg.framework/Versions/4/QtSvg48 \
	"${ImarisPath}/Contents/Frameworks/libbpfileio.7.6.dylib"
sudo install_name_tool -change @loader_path/QtSvg.framework/Versions/4/QtSvg \
	@loader_path/QtSvg.framework/Versions/4/QtSvg48 \
	"${ImarisPath}/Contents/Frameworks/libexpat.1.dylib"
sudo install_name_tool -change @loader_path/QtSvg.framework/Versions/4/QtSvg \
	@loader_path/QtSvg.framework/Versions/4/QtSvg48 \
	"${ImarisPath}/Contents/Frameworks/libfreeimage-3.15.0.dylib"
sudo install_name_tool -change @loader_path/QtSvg.framework/Versions/4/QtSvg \
	@loader_path/QtSvg.framework/Versions/4/QtSvg48 \
	"${ImarisPath}/Contents/Frameworks/libhdf5.7.dylib"
sudo install_name_tool -change @loader_path/QtSvg.framework/Versions/4/QtSvg \
	@loader_path/QtSvg.framework/Versions/4/QtSvg48 \
	"${ImarisPath}/Contents/Frameworks/libswscale.dylib"
###
sudo install_name_tool -change @loader_path/QtXml.framework/Versions/4/QtXml \
	@loader_path/QtXml.framework/Versions/4/QtXml48 \
	"${ImarisPath}/Contents/Frameworks/libCg.2.1.dylib"
sudo install_name_tool -change @loader_path/QtXml.framework/Versions/4/QtXml \
	@loader_path/QtXml.framework/Versions/4/QtXml48 \
	"${ImarisPath}/Contents/Frameworks/libCoin.60.dylib"
sudo install_name_tool -change @loader_path/QtXml.framework/Versions/4/QtXml \
	@loader_path/QtXml.framework/Versions/4/QtXml48 \
	"${ImarisPath}/Contents/Frameworks/libIce.34.dylib"
sudo install_name_tool -change @loader_path/QtXml.framework/Versions/4/QtXml \
	@loader_path/QtXml.framework/Versions/4/QtXml48 \
	"${ImarisPath}/Contents/Frameworks/libSBReadFile.dylib"
sudo install_name_tool -change @loader_path/QtXml.framework/Versions/4/QtXml \
	@loader_path/QtXml.framework/Versions/4/QtXml48 \
	"${ImarisPath}/Contents/Frameworks/libSoQt.20.dylib"
sudo install_name_tool -change @loader_path/QtXml.framework/Versions/4/QtXml \
	@loader_path/QtXml.framework/Versions/4/QtXml48 \
	"${ImarisPath}/Contents/Frameworks/libavcodec.dylib"
sudo install_name_tool -change @loader_path/QtXml.framework/Versions/4/QtXml \
	@loader_path/QtXml.framework/Versions/4/QtXml48 \
	"${ImarisPath}/Contents/Frameworks/libavformat.dylib"
sudo install_name_tool -change @loader_path/QtXml.framework/Versions/4/QtXml \
	@loader_path/QtXml.framework/Versions/4/QtXml48 \
	"${ImarisPath}/Contents/Frameworks/libavutil.dylib"
sudo install_name_tool -change @loader_path/QtXml.framework/Versions/4/QtXml \
	@loader_path/QtXml.framework/Versions/4/QtXml48 \
	"${ImarisPath}/Contents/Frameworks/libbpBatchIcePlugin.7.6.dylib"
sudo install_name_tool -change @loader_path/QtXml.framework/Versions/4/QtXml \
	@loader_path/QtXml.framework/Versions/4/QtXml48 \
	"${ImarisPath}/Contents/Frameworks/libbpImarisBatchAgent.7.6.dylib"
sudo install_name_tool -change @loader_path/QtXml.framework/Versions/4/QtXml \
	@loader_path/QtXml.framework/Versions/4/QtXml48 \
	"${ImarisPath}/Contents/Frameworks/libbpcore.7.6.dylib"
sudo install_name_tool -change @loader_path/QtXml.framework/Versions/4/QtXml \
	@loader_path/QtXml.framework/Versions/4/QtXml48 \
	"${ImarisPath}/Contents/Frameworks/libbpexception.7.6.dylib"
sudo install_name_tool -change @loader_path/QtXml.framework/Versions/4/QtXml \
	@loader_path/QtXml.framework/Versions/4/QtXml48 \
	"${ImarisPath}/Contents/Frameworks/libbpfileio.7.6.dylib"
sudo install_name_tool -change @loader_path/QtXml.framework/Versions/4/QtXml \
	@loader_path/QtXml.framework/Versions/4/QtXml48 \
	"${ImarisPath}/Contents/Frameworks/libexpat.1.dylib"
sudo install_name_tool -change @loader_path/QtXml.framework/Versions/4/QtXml \
	@loader_path/QtXml.framework/Versions/4/QtXml48 \
	"${ImarisPath}/Contents/Frameworks/libfreeimage-3.15.0.dylib"
sudo install_name_tool -change @loader_path/QtXml.framework/Versions/4/QtXml \
	@loader_path/QtXml.framework/Versions/4/QtXml48 \
	"${ImarisPath}/Contents/Frameworks/libhdf5.7.dylib"
sudo install_name_tool -change @loader_path/QtXml.framework/Versions/4/QtXml \
	@loader_path/QtXml.framework/Versions/4/QtXml48 \
	"${ImarisPath}/Contents/Frameworks/libswscale.dylib"
###
sudo install_name_tool -change @loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns \
	@loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns48 \
	"${ImarisPath}/Contents/Frameworks/libCg.2.1.dylib"
sudo install_name_tool -change @loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns \
	@loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns48 \
	"${ImarisPath}/Contents/Frameworks/libCoin.60.dylib"
sudo install_name_tool -change @loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns \
	@loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns48 \
	"${ImarisPath}/Contents/Frameworks/libIce.34.dylib"
sudo install_name_tool -change @loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns \
	@loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns48 \
	"${ImarisPath}/Contents/Frameworks/libSBReadFile.dylib"
sudo install_name_tool -change @loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns \
	@loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns48 \
	"${ImarisPath}/Contents/Frameworks/libSoQt.20.dylib"
sudo install_name_tool -change @loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns \
	@loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns48 \
	"${ImarisPath}/Contents/Frameworks/libavcodec.dylib"
sudo install_name_tool -change @loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns \
	@loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns48 \
	"${ImarisPath}/Contents/Frameworks/libavformat.dylib"
sudo install_name_tool -change @loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns \
	@loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns48 \
	"${ImarisPath}/Contents/Frameworks/libavutil.dylib"
sudo install_name_tool -change @loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns \
	@loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns48 \
	"${ImarisPath}/Contents/Frameworks/libbpBatchIcePlugin.7.6.dylib"
sudo install_name_tool -change @loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns \
	@loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns48 \
	"${ImarisPath}/Contents/Frameworks/libbpImarisBatchAgent.7.6.dylib"
sudo install_name_tool -change @loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns \
	@loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns48 \
	"${ImarisPath}/Contents/Frameworks/libbpcore.7.6.dylib"
sudo install_name_tool -change @loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns \
	@loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns48 \
	"${ImarisPath}/Contents/Frameworks/libbpexception.7.6.dylib"
sudo install_name_tool -change @loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns \
	@loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns48 \
	"${ImarisPath}/Contents/Frameworks/libbpfileio.7.6.dylib"
sudo install_name_tool -change @loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns \
	@loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns48 \
	"${ImarisPath}/Contents/Frameworks/libexpat.1.dylib"
sudo install_name_tool -change @loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns \
	@loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns48 \
	"${ImarisPath}/Contents/Frameworks/libfreeimage-3.15.0.dylib"
sudo install_name_tool -change @loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns \
	@loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns48 \
	"${ImarisPath}/Contents/Frameworks/libhdf5.7.dylib"
sudo install_name_tool -change @loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns \
	@loader_path/QtXmlPatterns.framework/Versions/4/QtXmlPatterns48 \
	"${ImarisPath}/Contents/Frameworks/libswscale.dylib"

# Change the references of other executables (in MacOS)
sudo install_name_tool -change \
	@executable_path/../Frameworks/QtCore.framework/Versions/4/QtCore \
	@executable_path/../Frameworks/QtCore.framework/Versions/4/QtCore48 \
	"${ImarisPath}/Contents/MacOS/Imaris"
sudo install_name_tool -change \
	@executable_path/../Frameworks/QtDeclarative.framework/Versions/4/QtDeclarative \
	@executable_path/../Frameworks/QtDeclarative.framework/Versions/4/QtDeclarative48 \
	"${ImarisPath}/Contents/MacOS/Imaris"
sudo install_name_tool -change \
	@executable_path/../Frameworks/QtGui.framework/Versions/4/QtGui \
	@executable_path/../Frameworks/QtGui.framework/Versions/4/QtGui48 \
	"${ImarisPath}/Contents/MacOS/Imaris"
sudo install_name_tool -change \
	@executable_path/../Frameworks/QtNetwork.framework/Versions/4/QtNetwork \
	@executable_path/../Frameworks/QtNetwork.framework/Versions/4/QtNetwork48 \
	"${ImarisPath}/Contents/MacOS/Imaris"
sudo install_name_tool -change \
	@executable_path/../Frameworks/QtOpenGL.framework/Versions/4/QtOpenGL \
	@executable_path/../Frameworks/QtOpenGL.framework/Versions/4/QtOpenGL48 \
	"${ImarisPath}/Contents/MacOS/Imaris"
sudo install_name_tool -change \
	@executable_path/../Frameworks/QtScript.framework/Versions/4/QtScript \
	@executable_path/../Frameworks/QtScript.framework/Versions/4/QtScript48 \
	"${ImarisPath}/Contents/MacOS/Imaris"
sudo install_name_tool -change \
	@executable_path/../Frameworks/QtSql.framework/Versions/4/QtSql \
	@executable_path/../Frameworks/QtSql.framework/Versions/4/QtSql48 \
	"${ImarisPath}/Contents/MacOS/Imaris"
sudo install_name_tool -change \
	@executable_path/../Frameworks/QtSvg.framework/Versions/4/QtSvg \
	@executable_path/../Frameworks/QtSvg.framework/Versions/4/QtSvg48 \
	"${ImarisPath}/Contents/MacOS/Imaris"
sudo install_name_tool -change \
	@executable_path/../Frameworks/QtXml.framework/Versions/4/QtXml \
	@executable_path/../Frameworks/QtXml.framework/Versions/4/QtXml48 \
	"${ImarisPath}/Contents/MacOS/Imaris"
sudo install_name_tool -change \
	@executable_path/../Frameworks/QtXmlPatterns.framework/Versions/4/QtXmlPatterns \
	@executable_path/../Frameworks/QtXmlPatterns.framework/Versions/4/QtXmlPatterns48 \
	"${ImarisPath}/Contents/MacOS/Imaris"
