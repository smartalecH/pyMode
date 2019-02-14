

class Line():
    def __init__(self, nLeft, nRight, X1,Y1,X2,Y2, *args, **kwargs):
        self.nLeft = nLeft
        self.nRight = nRight
        self.X1 = X1
        self.X2 = X2
        self.Y1 = Y1
        self.Y2 = Y2
    
    def writeContents(self):
            return "l ({:e},{:e}) ({:e},{:e}) {:e} {:e} {:e} {:e}\n".format(
            self.nLeft.real,self.nLeft.imag,self.nRight.real,self.nRight.imag,self.X1,self.Y1,self.X2,self.Y2
        )

class Bezier():
    def __init__(self, nLeft, nRight, X1,Y1,X2,Y2, *args, **kwargs):
        self.nLeft = nLeft
        self.nRight = nRight
        self.X1 = X1
        self.X2 = X2
        self.Y1 = Y1
        self.Y2 = Y2
    
    def writeContents(self):
            return "l ({:e},{:e}) ({:e},{:e}) {:e} {:e} {:e} {:e}\n".format(
            self.nLeft.real,self.nLeft.imag,self.nRight.real,self.nRight.imag,self.X1,self.Y1,self.X2,self.Y2
        )

class Rectangle():
    def __init__(self,centerX, centerY, width, height, core, cladding, rc=0, *args, **kwargs):
        self.centerX = centerX
        self.centerY = centerY
        self.width  = width
        self.height = height
        self.core = core
        self.cladding = cladding
        self.rc = rc

    def writeContents(self,wavelength):
        Lines = [
            [[self.centerX-self.width/2+self.rc,self.centerY+self.thickness/2],[self.centerX+self.width/2-self.rc,self.centerY+self.thickness/2]], # top line
            [[self.centerX-self.width/2+self.rc,self.centerY-self.thickness/2],[self.centerX+self.width/2-self.rc,self.centerY-self.thickness/2]], # bottom line
            [[self.centerX-self.width/2,self.centerY-self.thickness/2+self.rc],[self.centerX-self.width/2,self.centerY+self.thickness/2-self.rc]], # left line
            [[self.centerX+self.width/2,self.centerY-self.thickness/2+self.rc],[self.centerX+self.width/2,self.centerY+self.thickness/2-self.rc]]  # right line
        ]

        Corners = [
            [[self.centerX-self.width/2,self.centerY+self.thickness/2-self.rc],[self.centerX-self.width/2,self.centerY+self.thickness/2],[self.centerX-self.width/2+self.rc,self.centerY+self.thickness/2]], # top left
            [[self.centerX-self.width/2,self.centerY-self.thickness/2+self.rc],[self.centerX-self.width/2,self.centerY-self.thickness/2],[self.centerX-self.width/2+self.rc,self.centerY-self.thickness/2]], # bottom left
            [[self.centerX+self.width/2-self.rc,self.centerY+self.thickness/2],[self.centerX+self.width/2,self.centerY+self.thickness/2],[self.centerX+self.width/2,self.centerY+self.thickness/2-self.rc]], # top right
            [[self.centerX+self.width/2-self.rc,self.centerY-self.thickness/2],[self.centerX+self.width/2,self.centerY-self.thickness/2],[self.centerX+self.width/2,self.centerY-self.thickness/2+self.rc]]  # bottom right
        ]

        
        Indices = [
            [self.cladding,core],      # top line
            [core,cladding],      # bottom line
            [cladding,core],      # left line
            [core,cladding]       # right line
        ]

        fileContents = ""
        for k in range(len(Lines)):

            X1 = Lines[k][0][0]
            Y1 = Lines[k][0][1]
            X2 = Lines[k][1][0]
            Y2 = Lines[k][1][1]
            nLeft  = Indices[k][0]
            nRight = Indices[k][1]
            currentLine = "l ({:e},{:e}) ({:e},{:e}) {:e} {:e} {:e} {:e}\n".format(
                nLeft.real,nLeft.imag,nRight.real,nRight.imag,X1,Y1,X2,Y2
            )
            fileContents += (currentLine)
        if rc > 0:
            for k in range(len(Lines)):

                X1 = Corners[k][0][0]
                Y1 = Corners[k][0][1]
                X2 = Corners[k][1][0]
                Y2 = Corners[k][1][1]
                X3 = Corners[k][2][0]
                Y3 = Corners[k][2][1]
                nLeft  = Indices[k][0]
                nRight = Indices[k][1]
                currentLine = "b ({:e},{:e}) ({:e},{:e}) {:e} {:e} {:e} {:e} {:e} {:e}\n".format(
                    nLeft.real,nLeft.imag,nRight.real,nRight.imag,X1,Y1,X2,Y2,X3,Y3
                )
                fileContents += (currentLine)
        return fileContents

class Trapezoid():
    def __init__(self, *args, **kwargs):
        return super().__init__(*args, **kwargs)