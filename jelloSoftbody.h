#pragma once
#ifndef _JELLOSOFTBODY_H_
#define _JELLOSOFTBODY_H_

// == Internal Includes ==

#include "openGL-headers.h"
#include "pic.h"
#include "jello.h"

// == External Includes ==

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

class JelloSoftbody
{
private:
    // == Internal Temp Values ==

    point _vectSrcNgbr;
    point _relativeVelocity;
    double _vectSrcNgbrLength;
    double _elasticForce;

    // == Default Values ==

    const double _restLengthStructural = 1.0 / 7.0;
    const double _restLengthDiagonal2D = _restLengthStructural * sqrt(2.0);
    const double _restLengthDiagonal3D = sqrt(_restLengthStructural * _restLengthStructural + _restLengthDiagonal2D * _restLengthDiagonal2D);
    const double _restLengthBend = 2.0 / 7.0;

    // == Helper Functions ==

    void normalize(point& src)
    {
        double length;
        pLENGTH(src, length);
        src.x /= length;
        src.y /= length;
        src.z /= length;
    }

	// == Physics Functions ==

    void hookForce(const point& src, const point& ngbr, double kHook, point& force)
    {
        // Equation:
        // F = (-kH * (|L| - rL)) * (L/|L|)
        pMAKE(0.0, 0.0, 0.0, force);

        // un-normalized vector
        pDIFFERENCE(src, ngbr, _vectSrcNgbr); 
        pLENGTH(_vectSrcNgbr, _vectSrcNgbrLength); 
        _elasticForce = (-kHook) * (_vectSrcNgbrLength - _restLengthStructural); 

        // normalized vector
        normalize(_vectSrcNgbr);
        pMULTIPLY(_vectSrcNgbr, _elasticForce, force);
    }

    void hookForceDiagonal2D(const point& src, const point& ngbr, double kHook, point& force)
    {
        // Equation:
        // F = (-kH * (|L| - rL)) * (L/|L|)
        pMAKE(0.0, 0.0, 0.0, force);

        // un-normalized vector
        pDIFFERENCE(src, ngbr, _vectSrcNgbr);
        pLENGTH(_vectSrcNgbr, _vectSrcNgbrLength);
        _elasticForce = (-kHook) * (_vectSrcNgbrLength - _restLengthDiagonal2D);

        // normalized vector
        normalize(_vectSrcNgbr);
        pMULTIPLY(_vectSrcNgbr, _elasticForce, force);
    }

    void hookForceDiagonal3D(const point& src, const point& ngbr, double kHook, point& force)
    {
        // Equation:
        // F = (-kH * (|L| - rL)) * (L/|L|)
        pMAKE(0.0, 0.0, 0.0, force);

        // un-normalized vector
        pDIFFERENCE(src, ngbr, _vectSrcNgbr);
        pLENGTH(_vectSrcNgbr, _vectSrcNgbrLength);
        _elasticForce = (-kHook) * (_vectSrcNgbrLength - _restLengthDiagonal3D);

        // normalized vector
        normalize(_vectSrcNgbr);
        pMULTIPLY(_vectSrcNgbr, _elasticForce, force);
    }

    void dampingForce(const point& src, const point& ngbr, const point& srcVelocity, const point& ngbrVelocity, double kDamp, point& force)
    {
        // Equation:
        // F = (-kD * ((vA-vB) dot L) / |L|) * (L/|L|)
        pMAKE(0.0, 0.0, 0.0, force);

        // un-normalized vector
        pDIFFERENCE(src, ngbr, _vectSrcNgbr);    
        pDIFFERENCE(srcVelocity, ngbrVelocity, _relativeVelocity);
        pLENGTH(_vectSrcNgbr, _vectSrcNgbrLength);
        pDOT(_relativeVelocity, _vectSrcNgbr, _elasticForce);

        // normalized vector
        normalize(_vectSrcNgbr);
        _elasticForce *= (-kDamp) / _vectSrcNgbrLength; // (-kD * ((srcV - ngbrV) dot length) / |length|)
        pMULTIPLY(_vectSrcNgbr, _elasticForce, force);
    }

    void bendHookForce(const point& src, const point& ngbr, double kHook, point& force)
    {
        // Equation:
        // F = (-kH * (|L| - rL)) * (L/|L|)
        pMAKE(0.0, 0.0, 0.0, force);

        // un-normalized vector
        pDIFFERENCE(src, ngbr, _vectSrcNgbr);
        pLENGTH(_vectSrcNgbr, _vectSrcNgbrLength);
        _elasticForce = (-kHook) * (_vectSrcNgbrLength - _restLengthBend);

        // normalized vector
        normalize(_vectSrcNgbr);
        pMULTIPLY(_vectSrcNgbr, _elasticForce, force);
    }

    void hookForceCollision(const point& src, const point& ngbr, double kHook, point& normal, point& force)
    {
        // force clear and initalize
        pMAKE(0.0, 0.0, 0.0, force);
        // Calculate vector pointing from b to a
        pDIFFERENCE(src, ngbr, _vectSrcNgbr);
        pMULTIPLY(_vectSrcNgbr, -kHook, _vectSrcNgbr);
        pDOT(_vectSrcNgbr, normal, _elasticForce);
        pMULTIPLY(normal, _elasticForce, force);
    }

    void dampingForceCollision(const point& src, const point& ngbr, const point& velocitySrc, const point& velocityNgbr, double kDamp, point& n, point& force)
    {
        // force clear and initalize
        pMAKE(0.0, 0.0, 0.0, force);
        double _cosF;

        // un-normalized vector
        pDIFFERENCE(src, ngbr, _vectSrcNgbr);
        _vectSrcNgbrLength = sqrt(_vectSrcNgbr.x * _vectSrcNgbr.x + _vectSrcNgbr.y * _vectSrcNgbr.y + _vectSrcNgbr.z * _vectSrcNgbr.z);
        pDIFFERENCE(velocitySrc, velocityNgbr, _relativeVelocity);
        _elasticForce = (-kDamp) * ((_relativeVelocity.x * _vectSrcNgbr.x + _relativeVelocity.y * _vectSrcNgbr.y + _relativeVelocity.z * _vectSrcNgbr.z) / _vectSrcNgbrLength);
        
        // normalized vector
        normalize(_vectSrcNgbr);
        pMULTIPLY(_vectSrcNgbr, _elasticForce, force);
        pDOT(force, n, _cosF);
        pMULTIPLY(n, _cosF, force);
    }

public:
	// == Logic Functions == 

    /// <summary>
    /// Calculates the structural forces for mass located at jello->p[i][j][k]
    /// </summary>
    /// <param name="jello">: refrence object </param>
    /// <param name="i">: index </param>
    /// <param name="j">: index</param>
    /// <param name="k">: index</param>
    /// <param name="a">: acceleration value</param>
    void structuralForce(world* jello, int i, int j, int k, point& a)
    {
        point _force;

        if (i - 1 >= 0) 
        {
            hookForce(jello->p[i][j][k], jello->p[i - 1][j][k], jello->kElastic, _force);
            pSUM(_force, a, a);
            dampingForce(jello->p[i][j][k], jello->p[i - 1][j][k], jello->v[i][j][k], jello->v[i - 1][j][k], jello->dElastic, _force);
            pSUM(_force, a, a);
        }

        if (i + 1 <= 7)
        {
            hookForce(jello->p[i][j][k], jello->p[i + 1][j][k], jello->kElastic, _force);
            pSUM(_force, a, a);
            dampingForce(jello->p[i][j][k], jello->p[i + 1][j][k], jello->v[i][j][k], jello->v[i + 1][j][k], jello->dElastic, _force);
            pSUM(_force, a, a);
        }

        if (j - 1 >= 0)
        {
            hookForce(jello->p[i][j][k], jello->p[i][j - 1][k], jello->kElastic, _force);
            pSUM(_force, a, a);
            dampingForce(jello->p[i][j][k], jello->p[i][j - 1][k], jello->v[i][j][k], jello->v[i][j - 1][k], jello->dElastic, _force);
            pSUM(_force, a, a);
        }

        if (j + 1 <= 7)
        {
            hookForce(jello->p[i][j][k], jello->p[i][j + 1][k], jello->kElastic, _force);
            pSUM(_force, a, a);
            dampingForce(jello->p[i][j][k], jello->p[i][j + 1][k], jello->v[i][j][k], jello->v[i][j + 1][k], jello->dElastic, _force);
            pSUM(_force, a, a);
        }

        if (k - 1 >= 0)
        {
            hookForce(jello->p[i][j][k], jello->p[i][j][k - 1], jello->kElastic, _force);
            pSUM(_force, a, a);
            dampingForce(jello->p[i][j][k], jello->p[i][j][k - 1], jello->v[i][j][k], jello->v[i][j][k - 1], jello->dElastic, _force);
            pSUM(_force, a, a);
        }

        if (k + 1 <= 7)
        {
            hookForce(jello->p[i][j][k], jello->p[i][j][k + 1], jello->kElastic, _force);
            pSUM(_force, a, a);
            dampingForce(jello->p[i][j][k], jello->p[i][j][k + 1], jello->v[i][j][k], jello->v[i][j][k + 1], jello->dElastic, _force);
            pSUM(_force, a, a);
        }
    }

    /// <summary>
    /// Calculates the shear forces for mass located at jello->p[i][j][k]
    /// </summary>
    /// <param name="jello">: refrence object </param>
    /// <param name="i">: index </param>
    /// <param name="j">: index</param>
    /// <param name="k">: index</param>
    /// <param name="a">: acceleration value</param>
    void shearForce(world* jello, int i, int j, int k, point& a)
    {
        point _force;

        // 2D-3D space logic
        if (i > 0)
        {
            // 2d space
            if (k > 0) 
            {
                hookForceDiagonal2D(jello->p[i][j][k], jello->p[i - 1][j][k - 1], jello->kElastic, _force);
                pSUM(_force, a, a);
                dampingForce(jello->p[i][j][k], jello->p[i - 1][j][k - 1], jello->v[i][j][k], jello->v[i - 1][j][k - 1], jello->dElastic, _force);
                pSUM(_force, a, a);
            }

            if (k < 7)
            {
                hookForceDiagonal2D(jello->p[i][j][k], jello->p[i - 1][j][k + 1], jello->kElastic, _force);
                pSUM(_force, a, a);
                dampingForce(jello->p[i][j][k], jello->p[i - 1][j][k + 1], jello->v[i][j][k], jello->v[i - 1][j][k + 1], jello->dElastic, _force);
                pSUM(_force, a, a);
            }

            if (j > 0)
            {
                hookForceDiagonal2D(jello->p[i][j][k], jello->p[i - 1][j - 1][k], jello->kElastic, _force);
                pSUM(_force, a, a);
                dampingForce(jello->p[i][j][k], jello->p[i - 1][j - 1][k], jello->v[i][j][k], jello->v[i - 1][j - 1][k], jello->dElastic, _force);
                pSUM(_force, a, a);
            }

            if (j < 7) 
            {
                hookForceDiagonal2D(jello->p[i][j][k], jello->p[i - 1][j + 1][k], jello->kElastic, _force);
                pSUM(_force, a, a);
                dampingForce(jello->p[i][j][k], jello->p[i - 1][j + 1][k], jello->v[i][j][k], jello->v[i - 1][j + 1][k], jello->dElastic, _force);
                pSUM(_force, a, a);
            }
            // 3d space
            if (j > 0 && k > 0) 
            {
                hookForceDiagonal3D(jello->p[i][j][k], jello->p[i - 1][j - 1][k - 1], jello->kElastic, _force);
                pSUM(_force, a, a);
                dampingForce(jello->p[i][j][k], jello->p[i - 1][j - 1][k - 1], jello->v[i][j][k], jello->v[i - 1][j - 1][k - 1], jello->dElastic, _force);
                pSUM(_force, a, a);
            }

            if (j > 0 && k < 7) 
            {
                hookForceDiagonal3D(jello->p[i][j][k], jello->p[i - 1][j - 1][k + 1], jello->kElastic, _force);
                pSUM(_force, a, a);
                dampingForce(jello->p[i][j][k], jello->p[i - 1][j - 1][k + 1], jello->v[i][j][k], jello->v[i - 1][j - 1][k + 1], jello->dElastic, _force);
                pSUM(_force, a, a);
            }

            if (j < 7 && k > 0)
            {
                hookForceDiagonal3D(jello->p[i][j][k], jello->p[i - 1][j + 1][k - 1], jello->kElastic, _force);
                pSUM(_force, a, a);
                dampingForce(jello->p[i][j][k], jello->p[i - 1][j + 1][k - 1], jello->v[i][j][k], jello->v[i - 1][j + 1][k - 1], jello->dElastic, _force);
                pSUM(_force, a, a);
            }

            if (j < 7 && k < 7) 
            {
                hookForceDiagonal3D(jello->p[i][j][k], jello->p[i - 1][j + 1][k + 1], jello->kElastic, _force);
                pSUM(_force, a, a);
                dampingForce(jello->p[i][j][k], jello->p[i - 1][j + 1][k + 1], jello->v[i][j][k], jello->v[i - 1][j + 1][k + 1], jello->dElastic, _force);
                pSUM(_force, a, a);
            }
        }

        if (i < 7)
        {
            // 2d space
            if (k > 0) {
                hookForceDiagonal2D(jello->p[i][j][k], jello->p[i + 1][j][k - 1], jello->kElastic, _force);
                pSUM(_force, a, a);
                dampingForce(jello->p[i][j][k], jello->p[i + 1][j][k - 1], jello->v[i][j][k], jello->v[i + 1][j][k - 1], jello->dElastic, _force);
                pSUM(_force, a, a);
            }
            
            if (k < 7)
            {
                hookForceDiagonal2D(jello->p[i][j][k], jello->p[i + 1][j][k + 1], jello->kElastic, _force);
                pSUM(_force, a, a);
                dampingForce(jello->p[i][j][k], jello->p[i + 1][j][k + 1], jello->v[i][j][k], jello->v[i + 1][j][k + 1], jello->dElastic, _force);
                pSUM(_force, a, a);
            }

            if (j > 0)
            {
                hookForceDiagonal2D(jello->p[i][j][k], jello->p[i + 1][j - 1][k], jello->kElastic, _force);
                pSUM(_force, a, a);
                dampingForce(jello->p[i][j][k], jello->p[i + 1][j - 1][k], jello->v[i][j][k], jello->v[i + 1][j - 1][k], jello->dElastic, _force);
                pSUM(_force, a, a);
            }

            if (j < 7)
            {
                hookForceDiagonal2D(jello->p[i][j][k], jello->p[i + 1][j + 1][k], jello->kElastic, _force);
                pSUM(_force, a, a);
                dampingForce(jello->p[i][j][k], jello->p[i + 1][j + 1][k], jello->v[i][j][k], jello->v[i + 1][j + 1][k], jello->dElastic, _force);
                pSUM(_force, a, a);
            }
            // 3d cube
            if (j > 0 && k > 0)
            {
                hookForceDiagonal3D(jello->p[i][j][k], jello->p[i + 1][j - 1][k - 1], jello->kElastic, _force);
                pSUM(_force, a, a);
                dampingForce(jello->p[i][j][k], jello->p[i + 1][j - 1][k - 1], jello->v[i][j][k], jello->v[i + 1][j - 1][k - 1], jello->dElastic, _force);
                pSUM(_force, a, a);
            }

            if (j > 0 && k < 7)
            {
                hookForceDiagonal3D(jello->p[i][j][k], jello->p[i + 1][j - 1][k + 1], jello->kElastic, _force);
                pSUM(_force, a, a);
                dampingForce(jello->p[i][j][k], jello->p[i + 1][j - 1][k + 1], jello->v[i][j][k], jello->v[i + 1][j - 1][k + 1], jello->dElastic, _force);
                pSUM(_force, a, a);
            }

            if (j < 7 && k > 0)
            {
                hookForceDiagonal3D(jello->p[i][j][k], jello->p[i + 1][j + 1][k - 1], jello->kElastic, _force);
                pSUM(_force, a, a);
                dampingForce(jello->p[i][j][k], jello->p[i + 1][j + 1][k - 1], jello->v[i][j][k], jello->v[i + 1][j + 1][k - 1], jello->dElastic, _force);
                pSUM(_force, a, a);
            }

            if (j < 7 && k < 7)
            {
                hookForceDiagonal3D(jello->p[i][j][k], jello->p[i + 1][j + 1][k + 1], jello->kElastic, _force);
                pSUM(_force, a, a);
                dampingForce(jello->p[i][j][k], jello->p[i + 1][j + 1][k + 1], jello->v[i][j][k], jello->v[i + 1][j + 1][k + 1], jello->dElastic, _force);
                pSUM(_force, a, a);
            }
        }

        // 2D space logic
        if (j > 0)
        {
            // 2d space
            if (k > 0)
            {
                hookForceDiagonal2D(jello->p[i][j][k], jello->p[i][j - 1][k - 1], jello->kElastic, _force);
                pSUM(_force, a, a);
                dampingForce(jello->p[i][j][k], jello->p[i][j - 1][k - 1], jello->v[i][j][k], jello->v[i][j - 1][k - 1], jello->dElastic, _force);
                pSUM(_force, a, a);
            }

            if (k < 7)
            {
                hookForceDiagonal2D(jello->p[i][j][k], jello->p[i][j - 1][k + 1], jello->kElastic, _force);
                pSUM(_force, a, a);
                dampingForce(jello->p[i][j][k], jello->p[i][j - 1][k + 1], jello->v[i][j][k], jello->v[i][j - 1][k + 1], jello->dElastic, _force);
                pSUM(_force, a, a);
            }
        }

        if (j < 7) {
            // on 2D surface
            if (k > 0)
            {
                hookForceDiagonal2D(jello->p[i][j][k], jello->p[i][j + 1][k - 1], jello->kElastic, _force);
                pSUM(_force, a, a);
                dampingForce(jello->p[i][j][k], jello->p[i][j + 1][k - 1], jello->v[i][j][k], jello->v[i][j + 1][k - 1], jello->dElastic, _force);
                pSUM(_force, a, a);
            }

            if (k < 7)
            {
                hookForceDiagonal2D(jello->p[i][j][k], jello->p[i][j + 1][k + 1], jello->kElastic, _force);
                pSUM(_force, a, a);
                dampingForce(jello->p[i][j][k], jello->p[i][j + 1][k + 1], jello->v[i][j][k], jello->v[i][j + 1][k + 1], jello->dElastic, _force);
                pSUM(_force, a, a);
            }
        }
    }

    /// <summary>
    /// Calculates the bend forces for mass located at jello->p[i][j][k]
    /// </summary>
    /// <param name="jello">: refrence object </param>
    /// <param name="i">: index </param>
    /// <param name="j">: index</param>
    /// <param name="k">: index</param>
    /// <param name="a">: acceleration value</param>
    void bendForce(world* jello, int i, int j, int k, point& a) 
    {
        point _force;

        if (i > 1)
        {
            bendHookForce(jello->p[i][j][k], jello->p[i - 2][j][k], jello->kElastic, _force);
            pSUM(_force, a, a);
            dampingForce(jello->p[i][j][k], jello->p[i - 2][j][k], jello->v[i][j][k], jello->v[i - 2][j][k], jello->dElastic, _force);
            pSUM(_force, a, a);
        }

        if (i < 6)
        {
            bendHookForce(jello->p[i][j][k], jello->p[i + 2][j][k], jello->kElastic, _force);
            pSUM(_force, a, a);
            dampingForce(jello->p[i][j][k], jello->p[i + 2][j][k], jello->v[i][j][k], jello->v[i + 2][j][k], jello->dElastic, _force);
            pSUM(_force, a, a);
        }

        if (j > 1)
        {
            bendHookForce(jello->p[i][j][k], jello->p[i][j - 2][k], jello->kElastic, _force);
            pSUM(_force, a, a);
            dampingForce(jello->p[i][j][k], jello->p[i][j - 2][k], jello->v[i][j][k], jello->v[i][j - 2][k], jello->dElastic, _force);
            pSUM(_force, a, a);
        }

        if (j < 6)
        {
            bendHookForce(jello->p[i][j][k], jello->p[i][j + 2][k], jello->kElastic, _force);
            pSUM(_force, a, a);
            dampingForce(jello->p[i][j][k], jello->p[i][j + 2][k], jello->v[i][j][k], jello->v[i][j + 2][k], jello->dElastic, _force);
            pSUM(_force, a, a);
        }

        if (k > 1)
        {
            bendHookForce(jello->p[i][j][k], jello->p[i][j][k - 2], jello->kElastic, _force);
            pSUM(_force, a, a);
            dampingForce(jello->p[i][j][k], jello->p[i][j][k - 2], jello->v[i][j][k], jello->v[i][j][k - 2], jello->dElastic, _force);
            pSUM(_force, a, a);
        }

        if (k < 6)
        {
            bendHookForce(jello->p[i][j][k], jello->p[i][j][k + 2], jello->kElastic, _force);
            pSUM(_force, a, a);
            dampingForce(jello->p[i][j][k], jello->p[i][j][k + 2], jello->v[i][j][k], jello->v[i][j][k + 2], jello->dElastic, _force);
            pSUM(_force, a, a);
        }
    }

    /// <summary>
    /// Calculates the collision forces for mass located at jello->p[i][j][k]
    /// </summary>
    /// <param name="jello">: refrence object </param>
    /// <param name="i">: index </param>
    /// <param name="j">: index</param>
    /// <param name="k">: index</param>
    /// <param name="a">: acceleration value</param>
    void collisionForce(world* jello, int i, int j, int k, point& a) 
    {
        if ((jello->p[i][j][k].x <= -2) || (jello->p[i][j][k].x >= 2) ||
            (jello->p[i][j][k].y <= -2) || (jello->p[i][j][k].y >= 2) ||
            (jello->p[i][j][k].z <= -2) || (jello->p[i][j][k].z >= 2)) {
            point _force;
            point _normal;
            pMAKE(0.0, 0.0, 0.0, _normal);

            // initialize the collision point on the wall
            point _wall;
            pMAKE(0.0, 0.0, 0.0, _wall);
            _wall = jello->p[i][j][k];

            point _vectWall;
            pMAKE(0.0, 0.0, 0.0, _vectWall);

            // generate normal and wall values based on index
            if (jello->p[i][j][k].x <= -2)
            {
                _normal.x = 1;
                _wall.x = -2;
            }

            if (jello->p[i][j][k].x >= 2)
            {
                _normal.x = -1;
                _wall.x = 2;
            }

            if (jello->p[i][j][k].y <= -2)
            {
                _normal.y = 1;
                _wall.y = -2;
            }

            if (jello->p[i][j][k].y >= 2)
            {
                _normal.y = -1;
                _wall.y = 2;
            }

            if (jello->p[i][j][k].z <= -2)
            {
                _normal.z = 1;
                _wall.z = -2;
            }

            if (jello->p[i][j][k].z >= 2)
            {
                _normal.z = -1;
                _wall.z = 2;
            }

            // pushback forces
            hookForceCollision(jello->p[i][j][k], _wall, jello->kCollision, _normal, _force);
            pSUM(_force, a, a);

            dampingForceCollision(jello->p[i][j][k], _wall, jello->v[i][j][k], _vectWall, jello->dCollision, _normal, _force);
            pSUM(_force, a, a);
        }

    }

    /// <summary>
    /// Calculates the force feild forces for mass located at jello->p[i][j][k]
    /// </summary>
    /// <param name="jello">: refrence object </param>
    /// <param name="i">: index </param>
    /// <param name="j">: index</param>
    /// <param name="k">: index</param>
    /// <param name="a">: acceleration value</param>
    void forceFeildForce(world* jello, int i, int j, int k, point& a) {
        // index for force feild
        int _x, _y, _z;

        // Points mapped by binary (x, y, z)
        point _f000, _f001;
        point _f010, _f011;
        point _f100, _f101;
        point _f110, _f111;
        double _posI, _posJ, _posZ;
        double _gridLength;
        point _externalForce;
        pMAKE(0.0, 0.0, 0.0, _externalForce);

        // generate index point
        _x = int((jello->p[i][j][k].x + 2) * (jello->resolution - 1) / 4);
        _y = int((jello->p[i][j][k].y + 2) * (jello->resolution - 1) / 4);
        _z = int((jello->p[i][j][k].z + 2) * (jello->resolution - 1) / 4);

        // clamp index
        if (_x == (jello->resolution - 1)) 
        {
            _x--;
        }
        if (_y == (jello->resolution - 1)) 
        {
            _y--;
        }
        if (_z == (jello->resolution - 1)) 
        {
            _z--;
        }
        // applies force if point is within bounding bounding box
        if (((_x >= 0) && (_x <= jello->resolution - 1)) && ((_y >= 0) && (_y <= jello->resolution - 1)) && ((_y >= 0) && (_y <= jello->resolution - 1))) 
        {
            // gather the forces at each point
            _f000 = jello->forceField[(_x * jello->resolution * jello->resolution + _y * jello->resolution + _z)];
            _f001 = jello->forceField[(_x * jello->resolution * jello->resolution + _y * jello->resolution + (_z + 1))];
            _f010 = jello->forceField[(_x * jello->resolution * jello->resolution + (_y + 1) * jello->resolution + _z)];
            _f011 = jello->forceField[(_x * jello->resolution * jello->resolution + (_y + 1) * jello->resolution + (_z + 1))];
            _f100 = jello->forceField[((_x + 1) * jello->resolution * jello->resolution + _y * jello->resolution + _z)];
            _f101 = jello->forceField[((_x + 1) * jello->resolution * jello->resolution + _y * jello->resolution + (_z + 1))];
            _f110 = jello->forceField[((_x + 1) * jello->resolution * jello->resolution + (_y + 1) * jello->resolution + _z)];
            _f111 = jello->forceField[((_x + 1) * jello->resolution * jello->resolution + (_y + 1) * jello->resolution + (_z + 1))];

            // 3D interpolation
            _gridLength = 1.0 * 4 / (jello->resolution - 1);
            _posI = (jello->p[i][j][k].x - (-2 + 1.0 * 4 * _x / (jello->resolution - 1))) / _gridLength;
            _posJ = (jello->p[i][j][k].y - (-2 + 1.0 * 4 * _y / (jello->resolution - 1))) / _gridLength;
            _posZ = (jello->p[i][j][k].z - (-2 + 1.0 * 4 * _z / (jello->resolution - 1))) / _gridLength;

            // weighted trilinear interpolation for force at i,j,k
            pMULTIPLY(_f000, (1 - _posI) * (1 - _posJ) * (1 - _posZ), _f000);
            pMULTIPLY(_f001, (1 - _posI) * (1 - _posJ) * _posZ, _f001);
            pMULTIPLY(_f010, (1 - _posI) * _posJ * (1 - _posZ), _f010);
            pMULTIPLY(_f011, (1 - _posI) * _posJ * _posZ, _f011);
            pMULTIPLY(_f100, _posI * (1 - _posJ) * (1 - _posZ), _f100);
            pMULTIPLY(_f101, _posI * (1 - _posJ) * _posZ, _f101);
            pMULTIPLY(_f110, _posI * _posJ * (1 - _posZ), _f110);
            pMULTIPLY(_f111, _posI * _posJ * _posZ, _f111);

            // combining forces
            pSUM(_externalForce, _f000, _externalForce);
            pSUM(_externalForce, _f001, _externalForce);
            pSUM(_externalForce, _f010, _externalForce);
            pSUM(_externalForce, _f011, _externalForce);
            pSUM(_externalForce, _f100, _externalForce);
            pSUM(_externalForce, _f101, _externalForce);
            pSUM(_externalForce, _f110, _externalForce);
            pSUM(_externalForce, _f111, _externalForce);

            // store in return value
            pMAKE(_externalForce.x, _externalForce.y, _externalForce.z, a)
        }
    }

    /// <summary>
    /// Checks if a point is behind a given cube
    /// </summary>
    /// <param name="jello">: refrence object </param>
    /// <param name="p">: point </param>
    /// <param name="noCollision">: collison</param>
    /// <param name="normal">: normal</param>
    void checkInclinedCollision(world* jello, point& p, bool& noCollision, point& normal) {

        if (jello->incPlanePresent == 0)
        {
            return;
        }

        double _point;
        
        // dot product
        _point = jello->a * p.x + jello->b * p.y + jello->c * p.z + jello->d;

        int _pointDirection = 0;
        if (_point > 0) {
            _pointDirection = 1;
        }
        else if (_point < 0) {
            _pointDirection = -1;
        }

        pMAKE(jello->a, jello->b, jello->c, normal);
        // dot product
        _point = jello->a * normal.x + jello->b * normal.y + jello->c * normal.z + jello->d;

        int _planeDirection = 0;
        if (_point > 0) {
            _planeDirection = 1;
        }
        else if (_point < 0) {
            _planeDirection = -1;
        }

        noCollision = (_planeDirection * _pointDirection) == 1;
    }

    /// <summary>
    /// Calculates a collsion force if a collision with the inclined plane is present
    /// </summary>
    /// <param name="jello">: refrence object </param>
    /// <param name="i">: index </param>
    /// <param name="j">: index </param>
    /// <param name="k">: index </param>
    /// <param name="noCollision">: if collision</param>
    /// <param name="collisionNormal">: normal of collision</param>
    /// <param name="noCollision">: acceleration</param>
    void inclinedCollisionForce(world* jello, int i, int j, int k, const bool& noCollision, const point& collisionNormal, point& a) {
        if (noCollision == true)
            return;
       
        // D = (ax + by + cz + d) / sqrt(a*a + b*b + c*c)
        
        // distance vector
        double distance;
        distance = (jello->a * jello->p[i][j][k].x + jello->b * jello->p[i][j][k].y + jello->c * jello->p[i][j][k].z + jello->d);
        distance = distance / sqrt(jello->a * jello->a + jello->b * jello->b + jello->c * jello->c);
        
        // normalized vector
        point _normalizedNormal;
        point _normal;
        double _normalLength;

        _normalLength = sqrt(collisionNormal.x * collisionNormal.x + collisionNormal.y * collisionNormal.y + collisionNormal.z * collisionNormal.z);
        _normalizedNormal.x = (collisionNormal.x / _normalLength);
        _normalizedNormal.y = (collisionNormal.y / _normalLength);
        _normalizedNormal.z = (collisionNormal.z / _normalLength);
        pMULTIPLY(_normalizedNormal, distance, _normalizedNormal);
        
        // Calculate the point on plane
        point _pointPlane;
        pDIFFERENCE(jello->p[i][j][k], _normalizedNormal, _pointPlane);
        
        // direction
        _normal.x = (collisionNormal.x / _normalLength);
        _normal.y = (collisionNormal.y / _normalLength);
        _normal.z = (collisionNormal.z / _normalLength);

        // calculate collision
        point _force;
        hookForceCollision(jello->p[i][j][k], _pointPlane, jello->kCollision, _normal, _force);
        pSUM(_force, a, a);

        point _velocityPlane;
        pMAKE(0.0, 0.0, 0.0, _velocityPlane);
        dampingForceCollision(jello->p[i][j][k], _pointPlane, jello->v[i][j][k], _velocityPlane, jello->dCollision, _normal, _force);
        pSUM(_force, a, a);
    }
};

#endif