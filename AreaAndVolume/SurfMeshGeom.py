#!/usr/bin/env python
## -*- coding: utf-8 -*-
import numpy as np
import math
import sys

class SurfMeshGeom(object):

    def __init__(self,nodes,faces):
        self._nodes = nodes
        self._faces = faces 

    def GetNumOfNodes(self):
        """
        number of nodes
        """
        return len(self._nodes)
    
    def GetNumOfFaces(self):
        """
        number of faces
        """
        return len(self._faces)
    
    def GetNumOfIsoNodes(self):
        """
        check the isolate node number
        """
        nnodes = self.GetNumOfNodes()
        niso = 0
        k = 0
        while k < nnodes:
            (i,j) = np.where(self._faces == k)
            if len(i) == 0:
                niso += 1
            k += 1
        return niso

    def GetGenus(self):
        """
        calculate genus that is the number of holes in a surface
        """
        nnodes = self.GetNumOfNodes()
        nfaces = self.GetNumOfFaces()
        genus = (nfaces/2 - nnodes + 2)/2
        return genus

    def GetNormOfFaces(self):
        """
        calculate the triangle unit normal vector
        """
        nfaces = self.GetNumOfFaces()
        fnorm = np.empty((nfaces,3),float)
        k = 0
        for face in self._faces:
            p = self._TriNode(face)
            fnorm[k] = np.cross(p[2]-p[0],p[1]-p[0])
            k += 1
        return fnorm

    def GetArea(self):
        """
        get mesh surface area
        """
        area = 0
        for face in self._faces:
            area += self._TriArea(self._TriNode(face))
        return abs(area)

    def GetVolume(self):
        """
        get mesh volume
        """
        volume = 0
        for face in self._faces:
            p = self._TriNode(face)
            a = np.cross(p[1] - p[0], p[2] - p[0])
            b = (p[0] + p[1] + p[2])/3
            volume += (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])/6
        return abs(volume)
 
    def GetRatio(self):
        """
        get ratio of all triangles 
        """
        ratios = np.empty(self.GetNumOfFaces(), np.float)
        i = 0
        for face in self._faces:
            p = self._TriNode(face)
            edge = self._TriEdge(p)
            ratios[i] = min(edge)/max(edge)
            i += 1
        return ratios
  
    def GetMinAngle(self):
        """
        get the minimum angle
        """
        minAngle = math.pi
        for face in self._faces:
            p = self._TriNode(face)
            angle = self._TriAngle(self._TriEdge(p))
            if min(angle) < minAngle:
                minAngle = min(angle)
        return minAngle*180/math.pi
 
    def GetMaxAngle(self):
        """
        get the maximum angle
        """
        maxAngle = 0
        for face in self._faces:
            p = self._TriNode(face)
            angle = self._TriAngle(self._TriEdge(p))   
            if max(angle) > maxAngle:
                maxAngle = max(angle)
        return maxAngle*180/math.pi
    
    def GetAngle(self):
        """
        get the angles of all triangles
        """
        angles = np.empty(self.GetNumOfFaces()*3, np.float)
        i = 0 
        for face in self._faces:
            p = self._TriNode(face)
            angle = self._TriAngle(self._TriEdge(p))
            angles[i] = angle[0]*180/math.pi
            angles[i+1] = angle[1]*180/math.pi
            angles[i+2] = angle[2]*180/math.pi
            i += 3
        return angles
    
    def GetGaussianCurvature(self):
        """
        get the gaussian curvature at all vertices
        formula K = (2*pi-sum(theta(i)))/(A/3)
        """
        nnodes = self.GetNumOfNodes()
        gcurv = np.empty(nnodes, np.float)
        k = 0
        while k < nnodes:
            A = 0.0
            theta = 0.0
            (i,j) = np.where(self._faces == k)
            ntri = len(i)
            ii = 0
            while ii < ntri:
                p = self._TriNode(self._faces[i[ii]])
                A += self._TriArea(p)/3
                angle = self._TriAngle(self._TriEdge(p))
                theta += angle[j[ii]]
                ii += 1
            gcurv[k] = (2*math.pi-theta)/A
            k += 1
        return gcurv
    
    def GetMeanCurvature(self):
        """
        get the average mean curvature at all vertices
        formula Hn = 1/(4*A)*sum((cot(alpha)+cot(beta))*(v-vj))
        """
        nnodes = self.GetNumOfNodes()
        mcurv = np.empty((nnodes,3), np.float)
        k = 0
        while k < nnodes:
            A = 0.0
            B = 0.0
            (i,j) = np.where(self._faces == k)
            ntri = len(i)
            arounds = np.array(self._faces[i])
            ii = 0
            v0 = arounds[0]
            p0 = self._TriNode(v0)
            Q0 = k
            Q1 = arounds[0,j[0]-1]
            while ii < ntri:
                (m,n) = np.where(arounds == Q1)
                if np.all(arounds[m[0]] == v0):
                    v1 = arounds[m[1]]
                else:
                    v1 = arounds[m[0]]
                p1 = self._TriNode(v1)
                jj = 0 
                while v0[jj] == Q0 or v0[jj] == Q1:
                    jj += 1
                kk =0
                while v1[kk] == Q0 or v1[kk] == Q1:
                    kk += 1
                angle0 = self._TriAngle(self._TriEdge(p0))
                angle1 = self._TriAngle(self._TriEdge(p1))
                B += (1/math.tan(angle0[jj])+1/math.tan(angle1[kk]))*(self._nodes[Q0]-self._nodes[Q1])
                A += self._TriArea(p1)
                p0 = p1
                v0 = v1
                Q1 = v1[kk]
                ii += 1
            mcurv[k] = B/(4*A)
            k += 1
        return mcurv
    
    def GetSolidAngle(self):
        """
        calculate the solid angle at all nodes
        formula omega = sum(dihedral angle) - (n-2)*pi
        """
        nnodes = self.GetNumOfNodes()
        SolidAngs = np.empty(nnodes, np.float)
        k = 0
        while k < nnodes:
            theta = 0.0
            (i,j) = np.where(self._faces == k)
            ntri = len(i)
            arounds = np.array(self._faces[i])
            ii = 0
            v0 = arounds[0]
            p0 = self._TriNode(v0)
            Q0 = k
            Q1 = arounds[0,j[0]-1]
            while ii < ntri:
                (m,n) = np.where(arounds == Q1)
                if np.all(arounds[m[0]] == v0):
                    v1 = arounds[m[1]]
                else:
                    v1 = arounds[m[0]]
                p1 = self._TriNode(v1)
                jj = 0 
                while v0[jj] == Q0 or v0[jj] == Q1:
                    jj += 1
                kk =0
                while v1[kk] == Q0 or v1[kk] == Q1:
                    kk += 1
                nv0 = self._TriNorm(p0)
                nv1 = self._TriNorm(p1)
                tv = p0[jj] - p1[kk]
                alpha = self._IncludedAngle(nv0,nv1)
                if sum(np.multiply(nv0,tv))*sum(np.multiply(nv1,tv)) > 0:
                    theta += alpha
                else:
                    theta += math.pi - alpha
                p0 = p1
                v0 = v1
                Q1 = v1[kk]
                ii += 1
            SolidAngs[k] = theta - (ntri - 2)*math.pi
            k += 1
        return SolidAngs       
    
  


#内部函数
    def _TriNode(self,face):
        p = np.array([self._nodes[face[0]], self._nodes[face[1]], self._nodes[face[2]]])
        return p
    
    def _TriArea(self, p):
        return math.sqrt(sum(np.cross(p[1] - p[0], p[2] - p[0])**2))/2

    def _TriEdge(self,p):
        edge = np.array([math.sqrt(sum((p[0]-p[1])**2)), math.sqrt(sum((p[1]-p[2])**2)), math.sqrt(sum((p[0]-p[2])**2))])
        return edge
    
    def _TriAngle(self,edge):
        angle = np.array([np.arccos((edge[0]**2+edge[2]**2-edge[1]**2)/(2*edge[0]*edge[2])),\
                          np.arccos((edge[0]**2+edge[1]**2-edge[2]**2)/(2*edge[0]*edge[1])),np.arccos((edge[1]**2+edge[2]**2-edge[0]**2)/(2*edge[1]*edge[2]))])
        return angle
    
    def _TriNorm(self,p):
        return np.cross(p[1]-p[0],p[2]-p[0])
    
    def _IncludedAngle(self,n1,n2):
        return np.arccos(sum(np.multiply(n1,n2))/(math.sqrt(sum(n1**2))*math.sqrt(sum(n2**2))))
   
    
FILE_PATH = sys.argv[1]

class TestSurfMesh(object):

    def __init__(self):
        self.f = open(FILE_PATH, 'r')
        lines = self.f.readlines()
        self.fileType = lines[0]
        self.nnode, self.nface, self.unknow = self.splitNodeLine(lines[1])
        self.nodes = np.array(self.readNodes(lines[2:2 + int(self.nnode)]))
        self.faces = np.array(self.readFaces(lines[2 + int(self.nnode):]))

    def splitNodeLine(self, line):
        a,b,c = line.strip().split()
        return float(a.strip()), float(b.strip()), float(c.strip())

    def splitFaceLine(self, line):
        a,b,c,d = line.strip().split()
        return int(b.strip()), int(c.strip()), int(d.strip())

    def readNodes(self, lines):
        nodes = []
        for line in lines:
            nodes.append(self.splitNodeLine(line))
        return nodes

    def readFaces(self, lines):
        faces = []
        for line in lines:
            faces.append(self.splitFaceLine(line))
        return faces


test = TestSurfMesh()
mytest = SurfMeshGeom(test.nodes,test.faces)
print 'Area:', mytest.GetArea()
print 'Volume:', mytest.GetVolume()