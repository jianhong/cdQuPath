// This script allows to load the celltype from a csv file output from detectCellType.R
//
// Script author: Jianhong Ou
// Script should work with QuPath v0.3.0 or newer
// This source code is licensed under the MIT license
// Warning: this script only works for under QuPath GUI

import java.io.BufferedReader
import java.io.FileReader
import qupath.lib.io.*
import qupath.lib.geom.Point2

def csvFilePath  = Dialogs.promptForFile(null)


// Read the CSV file
def csvData = []
try {
    def reader = new BufferedReader(new FileReader(csvFilePath))
    def line
    while ((line = reader.readLine()) != null) {
        csvData.add(line.split(','))
    }
    reader.close()
} catch (Exception e) {
    e.printStackTrace()
}

// Find the column indices for 'x', 'y', and 'CellType'
// default position id is the first two columns, celltype is the third column
// otherwise, find it by 'Object.ID' and 'CellType'
def xIndex = 0
def yIndex = 1
def classIndex = 2

def headerRow = csvData[0]
headerRow.eachWithIndex { value, index ->
    if (value.replaceAll("\"", "") == 'x') {
        xIndex = index
    } else if (value.replaceAll("\"", "") == 'y') {
        yIndex = index
    } else if (value.replaceAll("\"", "") == 'CellType') {
        classIndex = index
    }
}

println(xIndex)
println(yIndex)
println(classIndex)

def x_asy = []
def y_asy = []
def class_asy = []
csvData.eachWithIndex { row, rowIndex ->
    if (rowIndex > 0) {
        def x = row[xIndex].replaceAll("\"", "")
        def y = row[yIndex].replaceAll("\"", "")
        def objectClass = row[classIndex].replaceAll("\"", "")
        x_asy.push(Float.parseFloat(x))
        y_asy.push(Float.parseFloat(y))
        class_asy.push(objectClass)
    }
}
println 'Imported cell types'
// Print the 'Object ID' to 'Class' map

def pathClasses = getQuPath().getAvailablePathClasses()

def cells = getCellObjects()
/*
def hierarchy = QP.getCurrentHierarchy()
def pixelSize = QP.getCurrentImageData().getServer().getPixelCalibration().getAveragedPixelSizeMicrons()
*/
def getCellClass(x, y, x_asy, y_asy, class_asy, maxDist=1) {
    def p = new Point2(x, y)
    def minDist = maxDist
    def ct = 'undefined'
    for(var i=0; i<x_asy.size(); i++) {
        var d = p.distance(x_asy[i], y_asy[i])
        if(d<=minDist) {
           minDist = d
           ct = class_asy[i]
       }
    }
    return(ct)
}

def celltype = cells.collect {
    def newClassID = getCellClass(it.getROI().getCentroidX(), it.getROI().getCentroidY(),
                                  x_asy, y_asy, class_asy)
    def newClass = PathClass.fromString(newClassID)
    if (!pathClasses.contains(newClass)) {
        pathClasses.add(newClass)}
    PathObjects.createDetectionObject(it.getROI(), newClass, it.getMeasurementList())

}

addObjects(celltype)

fireHierarchyUpdate()

println "Done!"

