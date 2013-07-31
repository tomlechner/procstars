// Star Catalog Image Generator
// By Tom Lechner, tomlechner.com
// 2013
//
// This program creates star globes from the wonderful Tycho2 star catalog,
// and also the Principal Galaxy Catalog.
//
// You can run with or with the gui.
//
// Usage:
//  ./procstars  -c ../Tycho.dat  -w32768  -H 7 -m12 -S 10 -a 40 -g
//  ./procstars  -c ../Tycho.cat  -S 6 -m 12 -s 12 -w 32768
//  ./procstars  -c ../Tycho.cat  -S 6 -m 12 -s 6 -w 16384
//  ./procstars  -c ../Tycho.cat  -o stars.tif
//
//
//
// This is released with the MIT license.
// This was first inspired by a generator written by Nathan Bergey from 01-22-2010,
// which he released CC Attribution 3.0 as Processing code, which processed the 
// Tycho 2 catalog.
//
//
// Fyi:
// At 8192 px wide, 23 pixels make 1 degree. 1/3 of a pixel is one arc minute, 1/200 of a pixel is 1 arc second.
// At 32768,        91 pixels make 1 degree, 1.5 pixels make one arc minute. 
// A typical human can discerne down to about 1 arcminute.
// the ISS takes up about 30 arc seconds, and is max about -5.9 in magnitude
// The moon takes up about 30 arc minutes, and is about magnitude -12.74
// The sun is about magnitude -26.74, and takes up about 32 arc minutes.
// The dimmest objects observed so far in human visible wavelengths is about magnitude 36.
//
//
//
// The MIT License (MIT)
// Copyright (c) 2013 Tom Lechner
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
//
//


//Todo:
//  implement galaxy approximations with principal galaxy catalog's morphology field
//  project billboard stuff like nebulae, LMC and SMC, M31, space stations, etc
//  make halo scale with magnitude
//  find nebula catalog?
//  what are those huge red smears in the sky??
//  fake the big galaxies LMC+SMC, m31?





#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <stdint.h>

#include <lax/laxoptions.h>
#include <lax/rowframe.h>
#include <lax/laxutils.h>
#include <lax/buttondowninfo.h>
#include <lax/quickfileopen.h>
#include <lax/checkbox.h>

#include </usr/include/GraphicsMagick/Magick++.h>

#include <lax/laxoptions.h>
#include <lax/anxapp.h>
#include <lax/curvewindow.h>
//#include <lax/iconselector.h>
#include <lax/button.h>
#include <lax/lineinput.h>
#include <lax/strmanip.h>

#include "catalogs.h"


#define PI (M_PI)
#define TWO_PI (2*M_PI)

//map: scale from [omin,omax] to [nmin,nmax]
#define map(value, omin,omax, nmin,nmax) (nmin+(value-omin)/(omax-omin)*(nmax-nmin))
#define constrain(v, min,max) (v<min ? min : (v>max ? max : v))
#define radians(deg) (deg/180.*M_PI)
#define degrees(rad) (rad*180./M_PI)


#define DBG

using namespace std;
using namespace Magick;
using namespace Laxkit;
using namespace LaxFiles;


//options:


//stuff to make computations easier:
//double pixelwidth=1./width;





//--------------------------------- HaloWindow ---------------------------------
/*! \class HaloWindow
 * Control window for adjusting map of star color index to rgb color.
 */
class HaloWindow : public Laxkit::RowFrame
{
  public:
	RenderContext *context;
	LaxImage *halo;
	CurveWindow *blowout, *ramp;
	ButtonDownInfo buttondown;
	

	HaloWindow(RenderContext *cntxt);
	virtual ~HaloWindow();
	virtual int init();
	virtual void Refresh();
	virtual int Event(const EventData *e,const char *mes);
	virtual void Reset(int which=~0);
	virtual int LBDown(int x,int y,unsigned int state,int count,const LaxMouse *d);
    virtual int LBUp(int x,int y,unsigned int state,const LaxMouse *d);
	virtual int MouseMove(int x,int y,unsigned int state,const LaxMouse *d);

	virtual void UpdateHalo();
	virtual void send();
};

HaloWindow::HaloWindow(RenderContext *cntxt)
  : RowFrame(NULL,"Halo Window","Halo Window",ROWFRAME_ROWS|ANXWIN_DOUBLEBUFFER, 0,310,600,250,0, NULL,0,"update")
{
	context=cntxt;
	halo=create_new_image(150,150);
}

HaloWindow::~HaloWindow()
{
}

int HaloWindow::LBDown(int x,int y,unsigned int state,int count,const LaxMouse *d)
{
	int xx=wholelist.e[0]->x();
	int yy=wholelist.e[0]->y();
	int ww=wholelist.e[0]->w()-2*pad;
	int hh=wholelist.e[0]->h()-2*pad;

	if (x>xx && x<xx+ww && y>yy && y<yy+hh) {
		if (x<xx+ww/2) buttondown.down(d->id, LEFTBUTTON, x,y, -1);
		else buttondown.down(d->id, LEFTBUTTON, x,y, 1);
	}

	return 0;
}

int HaloWindow::LBUp(int x,int y,unsigned int state,const LaxMouse *d)
{
	buttondown.up(d->id, LEFTBUTTON);
	return 0;
}

int HaloWindow::MouseMove(int x,int y,unsigned int state,const LaxMouse *d)
{
	if (!buttondown.any()) return 0;
	int ox,oy;
	buttondown.move(d->id, x,y, &ox,&oy);

	int xx=wholelist.e[0]->x();
	int yy=wholelist.e[0]->y();
	int ww=wholelist.e[0]->w()-2*pad;
	int hh=wholelist.e[0]->h()-2*pad;
	ox-=xx+ww/2;
	oy-=yy+hh/2;
	x-=xx+ww/2;
	y-=yy+hh/2;
	
	double od=sqrt(ox*ox+oy*oy);
	double nd=sqrt(x*x+y*y);

	context->usehalo=(context->usehalo-1)*od/nd + 1;
	if (context->usehalo<1.0001) context->usehalo=1.0001; //otherwise occasional mysterious glitches
	UpdateHalo();
	send();
	needtodraw=1;
	return 0;
}

int HaloWindow::init()
{
	anXWindow *last=NULL;

	CurveWindow *win=NULL;
	
	//AddSpacer(75,25,50,50,  200,100,200,50, -1); //file info
	AddSpacer(200,100,200,50,  200,100,200,50, -1); //preview area


	last=win=ramp=new CurveWindow(NULL,"Ramp","Ramp",0, 5,5,500,500,0, last,object_id,"ramp",
						 			 "Halo Ramp", "radius",0,1,  "a",0,255);
	//radius->AddPoint(.7,255); // *** or whatever halostart is
	ramp->SetInfo(context->ramp);
	ramp->GetInfo()->curvetype=CurveInfo::Autosmooth;
	AddWin(win,1, 200,100,200,50,0,  200,100,200,50,0, -1);


	last=win=blowout=new CurveWindow(NULL,"Blowout","Blowout",0, 5,5,500,500,0, last,object_id,"blowout",
							"Halo Blowout", "opacity",0,1,  "color",0,1);
	blowout->SetInfo(context->blowout);
	blowout->GetInfo()->curvetype=CurveInfo::Autosmooth;
	AddWin(win,1, 200,100,200,50,0,  200,100,200,50,0, -1);

	UpdateHalo();


	last->CloseControlLoop();
	Sync(1);
	return 0;
}

/*! Reset to default colors.
 * which is &1 for red, &2 for green, &4 for blue
 */
void HaloWindow::Reset(int which)
{
	context->Reset(which);
	ramp->Needtodraw(1);
	blowout->Needtodraw(1);

	needtodraw=1;
}

void HaloWindow::Refresh()
{
	if (!needtodraw) return;
	if (arrangedstate!=1) Sync(0);

	int pad=10;
	int w=wholelist.e[0]->w()-2*pad;
	int h=wholelist.e[0]->h()-2*pad;

	if (h<w) w=h;
	foreground_color(0);
	fill_rectangle(this,pad,pad,w,h);
	image_out_skewed(halo, this, pad,pad, w,0, 0,w);

	SwapBuffers();
	needtodraw=0;
}

void HaloWindow::send()
{
    if (win_owner) {
        SimpleMessage *ev=new SimpleMessage;
        app->SendMessage(ev,win_owner,win_sendthis,object_id);
    }
}

void HaloWindow::UpdateHalo()
{
	ramp->GetInfo()->RefreshLookup(256, 0,255);
	blowout->GetInfo()->RefreshLookup(256, 0,255);
	
	unsigned char *data=halo->getImageBuffer();
	int w=halo->w();

	//CreateStockHalo(w, context->usehalo, data, "c", ramp->GetInfo(),blowout->GetInfo(),"halo.png");
	CreateStockHalo(w, context->usehalo, data, "c", ramp->GetInfo(),blowout->GetInfo(),NULL);
	halo->doneWithBuffer(data);

	needtodraw=1;
}

int HaloWindow::Event(const EventData *e,const char *mes)
{
	const SimpleMessage *m=dynamic_cast<const SimpleMessage*>(e);
	if (!m) return anXWindow::Event(e,mes);

	if (!strcmp(mes,"ramp") || !strcmp(mes,"blowout")) {
		UpdateHalo();
		send();
	}

	needtodraw=1;
	return 0;
}




//--------------------------------- IndexWindow ---------------------------------
/*! \class IndexWindow
 * Control window for adjusting map of star color index to rgb color.
 */
class IndexWindow : public Laxkit::RowFrame
{
  public:
	RenderContext *context;
	CurveWindow *rr,*gg,*bb;

	IndexWindow(RenderContext *cntxt);
	virtual ~IndexWindow();
	virtual int init();
	virtual void Refresh();
	virtual int Event(const EventData *e,const char *mes);
	virtual void Reset(int which=~0);
	virtual void UpdateIndex();
	virtual void send();
};

IndexWindow::IndexWindow(RenderContext *cntxt)
  : RowFrame(NULL,"Index Window","Index Window",ROWFRAME_ROWS, 0,0,600,250,0, NULL,0,"update")
{
	context=cntxt;
}

IndexWindow::~IndexWindow()
{
}

int IndexWindow::init()
{
	anXWindow *last=NULL;

	CurveWindow *win=NULL;
	
	AddSpacer(75,25,50,50,  200,150,200,50, -1);


	last=win=rr=new CurveWindow(NULL,"Red","Red",0, 5,5,500,500,0, last,object_id,"red",
						 			 "Red", "B-V",-0.5,2,  "r",0,255);
	rr->SetInfo(context->index_r);
	AddWin(win,1, 200,100,200,50,0,  200,150,200,50,0, -1);


	last=win=gg=new CurveWindow(NULL,"Green","Green",0, 5,5,500,500,0, last,object_id,"green",
							"Green", "B-V",-0.5,2,  "g",0,255);
	gg->SetInfo(context->index_g);
	AddWin(win,1, 200,100,200,50,0,  200,150,200,50,0, -1);


	last=win=bb=new CurveWindow(NULL,"Blue","Blue",0, 5,5,500,500,0, last,object_id,"blue",
							"Blue", "B-V",-0.5,2, "b",0,255);
	bb->SetInfo(context->index_b);
	AddWin(win,1, 200,100,200,50,0,  200,150,200,50,0, -1);

	Reset(~0);

	last->CloseControlLoop();
	Sync(1);
	return 0;
}

/*! Reset to default colors.
 * which is &1 for red, &2 for green, &4 for blue
 */
void IndexWindow::Reset(int which)
{
	context->Reset(which);
	rr->Needtodraw(1);
	gg->Needtodraw(1);
	bb->Needtodraw(1);

	UpdateIndex();

	needtodraw=1;
}

void IndexWindow::Refresh()
{
	if (!needtodraw) return;
	if (arrangedstate!=1) Sync(0);

	int pad=20;
	int w=wholelist.e[0]->w()-2*pad;
	int h=wholelist.e[0]->h()-2*pad;
	double pos, r,g,b;

	StarColor color;
	for (int c=0; c<h; c++) {
		pos=(double)c/(h-1)*(2.5)-.5;

		//------
		indexToRgb(context,pos,0, color);
		r=color.redf();
		g=color.greenf();
		b=color.bluef();
		//------
		//pos=(double)c/(h-1)*255;
		//r=rr->GetInfo()->lookup[(int)pos]/255.;
		//g=gg->GetInfo()->lookup[(int)pos]/255.;
		//b=bb->GetInfo()->lookup[(int)pos]/255.;
		//------
		//r=rr->f(pos)/255;
		//g=gg->f(pos)/255;
		//b=bb->f(pos)/255;

		foreground_color(r,g,b);
		draw_line(this, pad,pad+c, pad+w,pad+c);
	}

	needtodraw=0;
}

void IndexWindow::send()
{
    if (win_owner) {
        SimpleMessage *ev=new SimpleMessage;
        app->SendMessage(ev,win_owner,win_sendthis,object_id);
    }
}

void IndexWindow::UpdateIndex()
{
	rr->GetInfo()->RefreshLookup(256, 0,255);
	gg->GetInfo()->RefreshLookup(256, 0,255);
	bb->GetInfo()->RefreshLookup(256, 0,255);
}

int IndexWindow::Event(const EventData *e,const char *mes)
{
	const SimpleMessage *m=dynamic_cast<const SimpleMessage*>(e);
	if (!m) return anXWindow::Event(e,mes);

	if (!strcmp(mes,"red") || !strcmp(mes,"green") || !strcmp(mes,"blue")) {
		send();
	}

	needtodraw=1;
	return 0;
}


//--------------------------------- MagnitudeSelectWindow ---------------------------------
/*! \class MagnitudeSelectWindow
 */
class MagnitudeSelectWindow : public Laxkit::anXWindow
{
  public:
	enum MagThings {
		Usemax=1,
		Usemin,
		Userange,
		Max,
		Min,
		Cutoff
	};

	int firsttime;
	int pad;
	RenderContext *context;
	ButtonDownInfo buttondown;

	double max, min;
	double cutoff;
	double usemax, usemin;

	MagnitudeSelectWindow(RenderContext *rr, anXWindow *parent);
	virtual ~MagnitudeSelectWindow() {}
	//virtual int init();
	virtual void Refresh();

	virtual int scan(int x,int y);
	virtual int pos(double mag);
	virtual double posToMag(int p);
	virtual void send();

	virtual int LBDown(int x,int y,unsigned int state,int count,const LaxMouse *d);
    virtual int LBUp(int x,int y,unsigned int state,const LaxMouse *d);
	virtual int MouseMove(int x,int y,unsigned int state,const LaxMouse *d);
};

MagnitudeSelectWindow::MagnitudeSelectWindow(RenderContext *rr, anXWindow *parent)
  : anXWindow(parent,"Mag Window","Mag Window",ANXWIN_DOUBLEBUFFER, 500,0,600,900,0, NULL,parent->object_id,"magselect")
{
	context=rr;
	pad=5;

	max=-1;
	min=15;
	cutoff=context->bigthreshhold;
	usemax=context->maximum_magnitude;
	usemin=context->minimum_magnitude;

	if (max<min) {
		if (usemin>min) usemin=min;
		if (usemax<max) usemax=max;
	} else {
		if (usemin<min) usemin=min;
		if (usemax>max) usemax=max;
	}

	win_colors=app->color_panel;
    win_colors->inc_count();
}

int MagnitudeSelectWindow::scan(int x,int y)
{
	if (x<win_w/2-pad) {
		int mx=pos(usemax), mn=pos(usemin);
		if (y>mx-2*pad && y<mx+2*pad) return Usemax;
		if (y>mn-2*pad && y<mn+2*pad) return Usemin;
		if (mn>mx) { int t=mn; mn=mx; mx=t; }
		if (y>mn && y<mx) return Userange;
	} else if (x>win_w/2+pad) {
		return Cutoff;
	} else {
		int textheight=app->defaultlaxfont->textheight();
		if (y<pad+textheight+pad) return Max;
		if (y>win_h-pad-textheight-pad) return Min;
	}

	return 0;
}

void MagnitudeSelectWindow::Refresh()
{
	if (!needtodraw) return;
	needtodraw=0;

	clear_window(this);

	int textheight=app->defaultlaxfont->textheight();
	char scratch[100];
	int p;

	//draw line between max and min
	double v;
	for (int c=pad+textheight+pad; c<win_h-2*pad-textheight; c++) {
		v=(double)(c-(pad+textheight+pad))/(win_h-4*pad-2*textheight);
		if (max<min) v=1-v;
		foreground_color(v,v,v);
		draw_line(this, win_w/2-pad/2,c, win_w/2+pad/2,c);
	}

	//draw line between usemax and usemin
	p=win_w/2-3*pad;
	int p2=pos(usemax), p1=pos(usemin);
	if (p1>p2) { int t=p1; p1=p2; p2=t; }
	for (int c=p1; c<p2; c++) {
		v=(double)(c-(pad+textheight+pad))/(win_h-4*pad-2*textheight);
		if (usemax<usemin) v=1-v;
		foreground_color(v,v,v);
		draw_line(this, p-pad/3,c, p+pad/3,c);
	}

	foreground_color(win_colors->fg);

	 //max
	sprintf(scratch,"%.7g",max);
	textout(this, scratch,-1, win_w/2,pad, LAX_HCENTER|LAX_TOP);

	 //min
	sprintf(scratch,"%.7g",min);
	textout(this, scratch,-1, win_w/2,win_h-pad, LAX_HCENTER|LAX_BOTTOM);

	 //usemax
	p=pos(usemax);
	sprintf(scratch,"%.3g",usemax);
	textout(this, scratch,-1, win_w/2-2*pad-textheight,p, LAX_RIGHT|LAX_VCENTER);
	draw_line(this, win_w/2-pad-textheight,p, win_w/2+pad,p);

	 //usemin
	p=pos(usemin);
	sprintf(scratch,"%.3g",usemin);
	textout(this, scratch,-1, win_w/2-2*pad-textheight,p, LAX_RIGHT|LAX_VCENTER);
	draw_line(this, win_w/2-pad-textheight,p, win_w/2+pad,p);

	
	 //cutoff
	p=pos(cutoff);
	sprintf(scratch,"%.3g",cutoff);
	draw_line(this, win_w/2-2*pad,p, win_w/2+pad+textheight,p);
	textout(this, scratch,-1, win_w/2+pad+textheight+pad,p, LAX_LEFT|LAX_VCENTER);

	textout(this, "halos",-1, win_w/2+pad,p-pad, LAX_LEFT|LAX_BOTTOM);
	textout(this, "points",-1, win_w/2+pad,p+pad, LAX_LEFT|LAX_TOP);

	SwapBuffers();
}

int MagnitudeSelectWindow::pos(double mag)
{
	int textheight=app->defaultlaxfont->textheight();
	return win_h-(2*pad+textheight) - (mag-min)/(max-min)*(win_h-4*pad-2*textheight);
}

double MagnitudeSelectWindow::posToMag(int p)
{
	int textheight=app->defaultlaxfont->textheight();
	return (win_h-(2*pad+textheight)-p)*(max-min) / (win_h-4*pad-2*textheight);
}

int MagnitudeSelectWindow::LBDown(int x,int y,unsigned int state,int count,const LaxMouse *d)
{
	int p=scan(x,y);
	if (p!=0) buttondown.down(d->id, LEFTBUTTON, x,y, p);
	return 0;
}

int MagnitudeSelectWindow::LBUp(int x,int y,unsigned int state,const LaxMouse *d)
{
	int dragged=buttondown.up(d->id, LEFTBUTTON);

	if (!dragged) {
		 //edit number in a box
	}

	return 0;
}

int MagnitudeSelectWindow::MouseMove(int x,int y,unsigned int state,const LaxMouse *d)
{
	int p=scan(x,y);
	cerr <<" mag select over: "<<p<<endl;

	int ox,oy;
	buttondown.move(d->id, x,y, &ox,&oy);

	if (!buttondown.any()) return 0;

	buttondown.getextrainfo(d->id,LEFTBUTTON, &p);

	double oldp=posToMag(oy);
	double newp=posToMag(y);

	if (p==Cutoff) {
		cutoff+=newp-oldp;
		if (cutoff>min) cutoff=min;
		if (cutoff<max) cutoff=max;
		context->bigthreshhold=cutoff;
		needtodraw=1;

		send();
		return 0;

	} else if (p==Usemin) {
		usemin+=newp-oldp;
		if (usemin<usemax) usemin=usemax;
		if (usemin>min) usemin=min;
		context->minimum_magnitude=usemin;
		needtodraw=1;

		send();
		return 0;

	} else if (p==Usemax) {
		usemax+=newp-oldp;
		if (usemax<max) usemax=max;
		if (usemax>usemin) usemax=usemin;
		context->maximum_magnitude=usemax;
		needtodraw=1;

		send();
		return 0;
	}


	return 0;
}

void MagnitudeSelectWindow::send()
{
    if (win_owner) {
        SimpleMessage *ev=new SimpleMessage;
        app->SendMessage(ev,win_owner,win_sendthis,object_id);
    }
}


//--------------------------------- SizeWindow ---------------------------------
/*! \class SizeWindow
 * Control window for adjusting sizes of stars.
 */
class SizeWindow : public Laxkit::RowFrame
{
  public:
	RenderContext *context;

	MagnitudeSelectWindow *magselect;
	CurveWindow *bigscale,*pointopacity;
	SizeWindow(RenderContext *rr);
	virtual ~SizeWindow();
	virtual int init();
	virtual void Refresh();
	virtual int Event(const EventData *e,const char *mes);
	virtual void Reset(int which=~0);
	virtual void UpdateSize();
	virtual void send(const char *mes=NULL);
};

SizeWindow::SizeWindow(RenderContext *rr)
  : RowFrame(NULL,"Size Window","Size Window",ROWFRAME_ROWS, 0,0,600,250,0, NULL,0,"update")
{
	context=rr;
	bigscale=NULL;
	pointopacity=NULL;
	magselect=NULL;
}

SizeWindow::~SizeWindow()
{
}

int SizeWindow::init()
{
	anXWindow *last=NULL;

	CurveWindow *win=NULL;
	
	MagnitudeSelectWindow *magselect=new MagnitudeSelectWindow(context, this);
	AddWin(magselect,1, 125,25,50,50,0,  200,150,200,50,0, -1);



	last=win=bigscale=new CurveWindow(NULL,"Large","Large",0, 5,5,500,500,0, last,object_id,"large",
						 			 "Big star scale", "magnitude",context->bigthreshhold,-1, "px",1,50);
	bigscale->SetInfo(context->bigscale);
	AddWin(win,1, 200,100,200,50,0,  200,150,200,50,0, -1);


	last=win=pointopacity=new CurveWindow(NULL,"Points","Points",0, 5,5,500,500,0, last,object_id,"points",
							"Point Opacity", "mag",context->minimum_magnitude,context->bigthreshhold,  "a",0,1);
	pointopacity->SetInfo(context->pointopacity);
	AddWin(win,1, 200,100,200,50,0,  200,150,200,50,0, -1);


//	last=win=pointamp=new CurveWindow(NULL,"Point Amp","Point Amp",0, 5,5,500,500,0, last,object_id,"pointamp",
//							"Point Amp", "mag",context->minimum_magnitude,context->bigthreshhold, "a",0,1);
//	context->pointamp=pointamp->GetInfo();
//	win->GetInfo()->curvetype=CurveInfo::Autosmooth;
//	win->MovePoint(1, 0,0);
//	AddWin(win,1, 200,100,200,50,0,  200,150,200,50,0, -1);

	UpdateSize();

	last->CloseControlLoop();
	Sync(1);
	return 0;
}

/*! Reset to default colors.
 * which is &1 for red, &2 for green, &4 for blue
 */
void SizeWindow::Reset(int which)
{
	context->Reset(which);
	bigscale->Needtodraw(1);
	pointopacity->Needtodraw(1);
	needtodraw=1;
}

void SizeWindow::Refresh()
{
	if (!needtodraw) return;
	if (arrangedstate!=1) Sync(0);

	needtodraw=0;
}

void SizeWindow::send(const char *mes)
{
    if (win_owner) {
        SimpleMessage *ev=new SimpleMessage;
        app->SendMessage(ev,win_owner,mes?mes:win_sendthis,object_id);
    }
}

void SizeWindow::UpdateSize()
{
	bigscale    ->GetInfo()->RefreshLookup(256, 1,bigscale->GetInfo()->ymax);
	pointopacity->GetInfo()->RefreshLookup(256, 0,255);
}

int SizeWindow::Event(const EventData *e,const char *mes)
{
	const SimpleMessage *m=dynamic_cast<const SimpleMessage*>(e);
	if (!m) return anXWindow::Event(e,mes);

	if (!strcmp(mes,"large") || !strcmp(mes,"points")) {
		UpdateSize();
		send();
	}

	if (!strcmp(mes,"magselect")) {
		bigscale->GetInfo()->SetXBounds(context->bigthreshhold,context->maximum_magnitude);
		bigscale->Needtodraw(1);
		pointopacity->GetInfo()->SetXBounds(context->minimum_magnitude,context->bigthreshhold);
		pointopacity->Needtodraw(1);
		send();
	}

	if (!strcmp(mes,"large")) {
		send("large");
	}


	needtodraw=1;
	return 0;
}


//--------------------------------- CatalogWindow ---------------------------------
/*! \class CatalogWindow
 */
class CatalogWindow : public Laxkit::RowFrame
{
  public:
	RenderContext *context;

	//MagnitudeSelectWindow *magselect;
	//CurveWindow *bigscale,*pointopacity;

	CatalogWindow(RenderContext *rr);
	virtual ~CatalogWindow();
	virtual int init();
	virtual void Refresh();
	virtual int Event(const EventData *e,const char *mes);
	virtual void Reset(int which=~0);
	virtual void UpdateCatalog();
	virtual void send(const char *mes=NULL);
};

CatalogWindow::CatalogWindow(RenderContext *rr)
  : RowFrame(NULL,"Catalog Window","Catalog Window",ROWFRAME_COLUMNS, 0,0,600,250,0, NULL,0,"update")
{
	context=rr;
}

CatalogWindow::~CatalogWindow()
{
}

int CatalogWindow::init()
{
	anXWindow *last=NULL;

	CheckBox *box;
	last=box=new CheckBox(this, "galactic","galactic",CHECK_LEFT, 0,0,0,0,0, last,win_owner,"galactic","Galactic");
	box->State(context->galactic?LAX_ON:LAX_OFF);
	AddWin(box,1, box->win_w,0,50,50,0,  box->win_h,10,10,50,0, -1);

	last=box=new CheckBox(this, "transparent","transparent",CHECK_LEFT, 0,0,0,0,0, last,win_owner,"transparent","Transparent");
	box->State(context->transparent?LAX_ON:LAX_OFF);
	AddWin(box,1, box->win_w,0,50,50,0,  box->win_h,10,10,50,0, -1);

	LineInput *line;
	last=line=new LineInput(this, "width","width",LINP_INT, 0,0,0,0,0, last,win_owner,"width","Width:");
	line->SetText((int)context->width);
	line->GetLineEdit()->win_style|=LINEEDIT_SEND_FOCUS_OFF;
	AddWin(line,1, line->win_w,0,300,50,0, line->win_h,0,10,50,0, -1);

	AddNull();

	AddSpacer(1000,900,50,50,  200,150,200,50, -1);


//	last=win=magmap=new CurveWindow(NULL,"Mag Map","Mag Map",0, 5,5,500,500,0, last,object_id,"magmap",
//							"Magnitude Map", "mag",context->minimum_magnitude,context->bigthreshhold, "a",0,1);
//	context->magmap=magmap->GetInfo();
//	win->GetInfo()->curvetype=CurveInfo::Autosmooth;
//	win->MovePoint(1, 0,0);
//	AddWin(win,1, 200,100,200,50,0,  200,150,200,50,0, -1);

	//UpdateCatalog();

	last->CloseControlLoop();
	Sync(1);
	return 0;
}

/*! Reset to default colors.
 * which is &1 for red, &2 for green, &4 for blue
 */
void CatalogWindow::Reset(int which)
{
	if (which&1) {
	}
	needtodraw=1;
}

void CatalogWindow::Refresh()
{
	if (!needtodraw) return;
	if (arrangedstate!=1) Sync(0);

	int pad=5;
	int xx=wholelist.e[4]->x()+pad;
	int yy=wholelist.e[4]->y()+pad;
	//int ww=wholelist.e[4]->w();
	//int hh=wholelist.e[4]->h();

	textout(this, "Catalogs:",-1, xx,yy,LAX_LEFT|LAX_TOP);

	char scratch[200];
	int textheight=app->defaultlaxfont->textheight();
	Catalog *cat;
	int y=yy+textheight;
	for (int c=0; c<context->catalogs.n; c++) {
		cat=context->catalogs.e[c];
		//sprintf(scratch,"^ v x %s: %s",cat->name,lax_basename(cat->filename));
		sprintf(scratch,"%s: %s",cat->name,lax_basename(cat->filename));
		textout(this, scratch,-1, xx,y,LAX_LEFT|LAX_TOP);
		y+=textheight;
	}
	//textout(this, "New: PGC  Tycho  Random  Galaxy",-1, xx,y,LAX_LEFT|LAX_TOP);

	needtodraw=0;
}

void CatalogWindow::send(const char *mes)
{
    if (win_owner) {
        SimpleMessage *ev=new SimpleMessage;
        app->SendMessage(ev,win_owner,mes?mes:win_sendthis,object_id);
    }
}

void CatalogWindow::UpdateCatalog()
{
	//bigscale    ->GetInfo()->RefreshLookup(256, 1,bigscale->GetInfo()->ymax);
	//pointopacity->GetInfo()->RefreshLookup(256, 0,255);
}

int CatalogWindow::Event(const EventData *e,const char *mes)
{
	const SimpleMessage *m=dynamic_cast<const SimpleMessage*>(e);
	if (!m) return anXWindow::Event(e,mes);

	if (!strcmp(mes,"large") || !strcmp(mes,"points")) {
		UpdateCatalog();
		send();
	}


	needtodraw=1;
	return 0;
}




//--------------------------------- MainWindow ---------------------------------
/*! \class MainWindow
 */
class MainWindow : public Laxkit::RowFrame
{
  public:
	int firsttime;
	int zoomed;

	int numstars;
	RandomCatalog previewcatalog;
	LaxImage *preview;
	int needtoupdate;
	ButtonDownInfo buttondown;

	LineInput *saveto;
	HaloWindow *halowindow;
	CatalogWindow *catalogwindow;

	LaxImage *transparentmask;

	RenderContext *context;
	MainWindow(RenderContext *rr);
	virtual ~MainWindow();
	virtual int init();
	virtual void Refresh();
	virtual int Event(const EventData *e,const char *mes);
	//virtual void Reset(int which=~0);
	virtual void UpdatePreview();

	virtual int scan(int x,int y);
	virtual int LBDown(int x,int y,unsigned int state,int count,const LaxMouse *d);
    virtual int LBUp(int x,int y,unsigned int state,const LaxMouse *d);
	virtual int MouseMove(int x,int y,unsigned int state,const LaxMouse *d);
};

MainWindow::MainWindow(RenderContext *rr)
  : RowFrame(NULL,"Main Window","Main Window",ANXWIN_ESCAPABLE|ROWFRAME_ROWS|ANXWIN_DOUBLEBUFFER, 500,0,600,900,0, NULL,0,NULL),
	previewcatalog("preview",1000)
{
	numstars=1000;
	context=rr;
	preview=NULL;
	firsttime=2;
	needtoupdate=1;
	zoomed=0;

	transparentmask=NULL;
}

MainWindow::~MainWindow()
{
	transparentmask->dec_count();
}

int MainWindow::init()
{
	anXWindow *last=NULL;

	AddSpacer(win_w,0,1000,50,  win_w/2,0,500,50, -1); //area for preview image
	//AddSpacer(win_w/2,0,1000,50,  win_w/2,0,500,50, -1); //area for preview image
	//AddSpacer(200,10,2000,50,  10,10,200,50, -1); //area for preview image
	AddNull();


	//  Save To: ____ [...]    <<RENDER NOW>>    Save Settings To: ______[...] [load]
	last=saveto=new LineInput(this, "File","File",0, 0,0,0,0,0, NULL,object_id,"file","Render To:",context->filename);
	AddWin(saveto,1, saveto->win_w,0,300,50,0, saveto->win_h,0,10,50,0, -1);

	last=new QuickFileOpen(this, "Saveto","Saveto",0, 0,0,0,0,0, last,object_id,"saveto", FILES_SAVE, saveto);
	AddWin(last,1, last->win_w,0,0,50,0, last->win_h*.9,0,0,50,0, -1);


	AddSpacer(10,10,2000,50,  10,10,200,50, -1); //area for preview image

	Button *button;
	last=button=new Button(this, "Render","Render",0, 0,0,0,0,0, last,object_id,"render",0, " Render now! ");
	button->gap=button->win_w*.25;
	AddWin(button,1, button->win_w,0,300,50,0, button->win_h,0,10,50,0, -1);

	AddSpacer(10,10,2000,50,  10,10,200,50, -1); //area for preview image


	MessageBar *m=new MessageBar(this, "Save Settings","Save Settings",0, 0,0,0,0,0, "Project: ");
	AddWin(m,1, m->win_w,0,0,50,0, m->win_h,0,0,50,0, -1);


	last=new Button(this, "Save Settings","Save Settings",0, 0,0,0,0,0, last,object_id,"saveprojto",0,"Save");
	AddWin(last,1, last->win_w,0,300,50,0, last->win_h,0,10,50,0, -1);

	last=new Button(this, "Load Settings","Load Settings",0, 0,0,0,0,0, last,object_id,"loadproj",0,"Load");
	AddWin(last,1, last->win_w,0,300,50,0, last->win_h,0,10,50,0, -1);

//	last=saveprojto=new LineInput(this, "Save Settings","Save Settings",0, 0,0,0,0,0, last,object_id,"settings","Project:",context->projectfile);
//	AddWin(saveprojto,1, saveprojto->win_w,0,300,50,0, saveprojto->win_h,0,10,50,0, -1);
//
//	last=new QuickFileOpen(this, "Saveprojto","Saveprojto",0, 0,0,0,0,0, last,object_id,"saveprojto", FILES_SAVE, saveprojto);
//	AddWin(last,1, last->win_w,0,0,50,0, last->win_h*.9,0,0,50,0, -1);

	AddNull();


	 //-----editing windows:
	SizeWindow *swin=new SizeWindow(context);
	swin->SetOwner(this);
	AddWin(swin,1, 800,400,200,50,0,  200,150,200,50,0, -1);

	IndexWindow *iwin=new IndexWindow(context);
	iwin->SetOwner(this);
	AddWin(iwin,1, 800,400,200,50,0,  200,150,200,50,0, -1);

	halowindow=new HaloWindow(context);
	halowindow->SetOwner(this);
	AddWin(halowindow,1, 800,400,200,50,0,  200,150,200,50,0, -1);

	catalogwindow=new CatalogWindow(context);
	catalogwindow->SetOwner(this);
	AddWin(catalogwindow,1, 800,400,200,50,0,  200,150,200,50,0, -1);

	//CatalogWindow *catwin=new CatalogWindow();
	//AddWin(catwin,1, 800,400,200,50,0,  200,150,200,50,0, -1);



	last->CloseControlLoop();
	Sync(1);


	int bw=50;
	transparentmask=create_new_image(bw,bw);
	unsigned char *data=transparentmask->getImageBuffer();
	for (int i=0; i<bw*bw*4; i+=4) {
		data[i+0]=0;
		data[i+1]=0;
		data[i+2]=0;
		data[i+3]=200;
	}
	transparentmask->doneWithBuffer(data);


	return 0;
}

int MainWindow::Event(const EventData *e,const char *mes)
{
	const SimpleMessage *m=dynamic_cast<const SimpleMessage*>(e);
	if (!m) return anXWindow::Event(e,mes);

	if (!strcmp(mes,"galactic")) {
		if (m->info1==LAX_ON) context->galactic=1;
		else context->galactic=0;
		needtodraw=1;
		UpdatePreview();
		return 0;
	}
	
	if (!strcmp(mes,"transparent")) {
		if (m->info1==LAX_ON) context->transparent=1;
		else context->transparent=0;
		needtodraw=1;
		return 0;
	}

	if (!strcmp(mes,"width")) {
		context->width=strtol(m->str,NULL,10);
		context->height=context->width/2;
		needtodraw=1;
		return 0;
	}



	if (!strcmp(mes,"saveto")) {
		makestr(context->filename,m->str);
		saveto->SetText(m->str);
		saveto->GetLineEdit()->SetCurpos(-1);
		return 0;
	}

	if (!strcmp(mes,"loadproj")) {
	    app->rundialog(new FileDialog(NULL,"Load Project","Load Project", 0,
                                0,0,0,0,1, 
                                object_id,"nowloadproj",
                                FILES_OPEN_ONE,
                                context->projectfile
                               ));
		return 0;
	}

	if (!strcmp(mes,"nowloadproj")) {
		FILE *f=fopen(context->projectfile,"r");
		if (!f) {
			cerr << "Could not save to "<<context->projectfile<<"!"<<endl;
		} else {
			makestr(context->projectfile,m->str);
			context->dump_in(f,0,0,NULL,NULL);
			fclose(f);

			UpdatePreview();
		}
		return 0;
	}

	if (!strcmp(mes,"saveprojto")) {
	    app->rundialog(new FileDialog(NULL,"Save Project","Save Project", 0,
                                0,0,0,0,1, 
                                object_id,"nowsaveprojto",
                                FILES_SAVE,
                                context->projectfile
                               ));
		return 0;
	}

	if (!strcmp(mes,"nowsaveprojto")) {
		makestr(context->projectfile,m->str);
		FILE *f=fopen(context->projectfile,"w");
		if (!f) {
			cerr << "Could not save to "<<context->projectfile<<"!"<<endl;
		} else {
			context->dump_out(f,0,0,NULL);
			fclose(f);
		}
		return 0;
	}

	if (!strcmp(mes,"render")) {
		cout <<"Render..."<<endl;

		Render(context);
		return 0;
	}

	if (!strcmp(mes,"large")) {
		 //max star size changed, need to update halo sample
		//halowindow->***;
	}

	if (!strcmp(mes,"update")) {
		needtoupdate=1;
		needtodraw=1;
		return 0;
	}

	return anXWindow::Event(e,mes);
}

void MainWindow::UpdatePreview()
{
	cerr <<"Updating preview..."<<endl;

	if (context->halo && context->halowidth<context->maxstarsize*context->usehalo) {
		 //star size was changed, we need to make a different size halo
		cerr << "Resizing halo sample size.."<<endl;
		delete[] context->halo;
		context->halo=NULL;
	}
	if (!context->halo) {
		halowindow->UpdateHalo();
		if (!context->halo) {
			context->halowidth=2*context->maxstarsize*context->usehalo;
			context->halo=new unsigned char[context->halowidth*context->halowidth];
		}
	}
	CreateStockHalo(context->halowidth, context->usehalo, context->halo, "g",
					context->ramp,context->blowout,NULL);
					//context->ramp,context->blowout,"halo.png");

	int ww=wholelist.e[0]->w();
	int hh=wholelist.e[0]->h();

	if (preview && (ww!=preview->w() || hh!=preview->h())) {
		preview->dec_count();
		preview=NULL;
	}
	if (!preview) {
		preview=create_new_image(ww,hh);
	}

	unsigned char *data=preview->getImageBuffer();
	ww=preview->w();
	hh=preview->h();
	memset(data,0,ww*hh*4);
	for (int x=0; x<ww; x++) {
	  for (int y=0; y<hh; y++) {
		data[(y*ww+x)*4+3]=255;
	  }
	}

	previewcatalog.Render(context, data,ww,hh);
	preview->doneWithBuffer(data);
}

void MainWindow::Refresh()
{
	if (!needtodraw) return;
	if (arrangedstate!=1) Sync(0);

	if (firsttime>0) {
		firsttime--;
		return;
	}


	needtodraw=0;

	int xx=wholelist.e[0]->x();
	int yy=wholelist.e[0]->y();
	int ww=wholelist.e[0]->w();
	int hh=wholelist.e[0]->h();

	if (!preview || needtoupdate) {
		UpdatePreview();
		needtoupdate=0;
	}

	image_out_skewed(preview, this, xx,yy, ww,0, 0,hh);

	foreground_color(.5,.5,.5);
	if (!zoomed) {
		
		int min_asc_p=xx+context->min_asc/360*ww;
		int max_asc_p=xx+context->max_asc/360*ww;
		int min_dec_p=yy+context->min_dec/180*hh+hh/2;
		int max_dec_p=yy+context->max_dec/180*hh+hh/2;


		 //draw blackout
		if (context->min_asc>0) {
			image_out_skewed(transparentmask, this, 0,0, min_asc_p,0, 0,hh);
		}
		if (context->max_asc<360) {
			image_out_skewed(transparentmask, this, max_asc_p,0, ww-max_asc_p+1,0, 0,hh);
		}
		if (context->min_dec>-90) {
			image_out_skewed(transparentmask, this, 0,0, ww,0, 0,min_dec_p);
		}
		if (context->max_dec<90) {
			image_out_skewed(transparentmask, this, 0,max_dec_p, ww,0, 0,hh-max_dec_p);
		}

		 //draw lines
		if (context->min_asc>0) {
			draw_line(this, min_asc_p,yy, min_asc_p,yy+hh);
		}
		if (context->max_asc<360) {
			draw_line(this, max_asc_p,yy, max_asc_p,yy+hh);
		}
		if (context->min_dec>-90) {
			draw_line(this, xx,min_dec_p, xx+ww,min_dec_p);
		}
		if (context->max_dec<90) {
			draw_line(this, xx,max_dec_p, xx+ww,max_dec_p);
		}
	}

	SwapBuffers();
}

#define MINASC 1
#define MAXASC 2
#define MINDEC 3
#define MAXDEC 4
#define SWITCH 5
#define REFRESH 6
#define FEWERSTARS 7
#define MORESTARS  8

int MainWindow::scan(int x,int y)
{
	int range=20;

	int xx=wholelist.e[0]->x();
	int yy=wholelist.e[0]->y();
	int ww=wholelist.e[0]->w();
	int hh=wholelist.e[0]->h();
	int p;

	if (!zoomed) {
		p=xx+context->min_asc/360*ww;
		cerr <<"x:"<<x<<"  p:"<<p<<"  minasc:"<<context->min_asc<<endl;
		if (x<p+range && x>p-range) return MINASC;

		p=xx+context->max_asc/360*ww;
		cerr <<"x:"<<x<<"  p:"<<p<<"  maxasc:"<<context->max_asc<<endl;
		if (x<p+range && x>p-range) return MAXASC;

		p=yy+context->min_dec/180*hh+hh/2;
		if (y<p+range && y>p-range) return MINDEC;

		p=yy+context->max_dec/180*hh+hh/2;
		if (y<p+range && y>p-range) return MAXDEC;
	}

	return -1;
}

int MainWindow::LBDown(int x,int y,unsigned int state,int count,const LaxMouse *d)
{
	int found=scan(x,y);

	if (found>0) {
		buttondown.down(d->id,LEFTBUTTON, x,y, found);
	}


	return 0;
}

int MainWindow::LBUp(int x,int y,unsigned int state,const LaxMouse *d)
{
	int found=-1;
	buttondown.up(d->id,LEFTBUTTON, &found);

	if (found==REFRESH || found==MORESTARS || found==FEWERSTARS) {
		if (x<win_w/3) numstars*=.6;
		else if (x>win_w*2/3) numstars*=1.3;
		if (numstars<100) numstars=100;

		cerr <<"Regenerate "<<numstars<<"stars..."<<endl;


		previewcatalog.Repopulate(numstars,0);
		UpdatePreview();
		needtodraw=1;
		return 0;
	}

	return 0;
}

int MainWindow::MouseMove(int x,int y,unsigned int state,const LaxMouse *d)
{
	DBG int ff=scan(x,y);
	DBG cerr <<"over: "<<ff<<endl;


	if (!buttondown.any()) return 0;
	
	int ox,oy;
	int found=-1;
	buttondown.getextrainfo(d->id,LEFTBUTTON, &found);
	if (found<=0) return 0;
	buttondown.move(d->id, x,y, &ox,&oy);

	int dx=x-ox;
	int dy=y-oy;

	int ww=wholelist.e[0]->w();
	int hh=wholelist.e[0]->h();

	double ddx=double(dx)/ww*360;
	double ddy=double(dy)/hh*180;

	if (found==MINASC) {
		context->min_asc+=ddx;
		if (context->min_asc>context->max_asc) context->min_asc=context->max_asc;
		else if (context->min_asc<0) context->min_asc=0;

	} else if (found==MAXASC) {
		context->max_asc+=ddx;
		if (context->max_asc<context->min_asc) context->max_asc=context->min_asc;
		else if (context->max_asc>360) context->max_asc=360;

	} else if (found==MINDEC) {
		context->min_dec+=ddy;
		if (context->min_dec>context->max_dec) context->min_dec=context->max_dec;
		else if (context->min_dec<-90) context->min_dec=-90;

	} else if (found==MAXDEC) {
		context->max_dec+=ddy;
		if (context->max_dec<context->min_dec) context->max_dec=context->min_dec;
		else if (context->max_dec>90) context->max_dec=90;
	}

	needtodraw=1;
	return 0;
}

//--------------------------------- main() --------------------------------

LaxOptions options;

void InitOptions()
{
	options.HelpHeader( "Star Catalog Image Generator\n"
						"By Tom Lechner, 2013\n"
						"Version 0.000000001");
    options.UsageLine("procstars [options]");

	options.Add("width",        'w', 1, "Width of equirectangular texture image",                      0, "8192");
	options.Add("output-file",  'o', 1, "Filename for generated texture image (will save as tif)",     0, "file");
	options.Add("catalog-file", 'c', 1, "File containing the Tycho 2 Star Catalog",                    0, "catalog");
	options.Add("pgc-file",     'p', 1, "File containing the Principal Galaxy Catalog",                0, "catalog");
	options.Add("magnitude",    'm', 1, "Minimum brightness. stars in catalog are 0 (bright) "
			                            "to 15 (very dim). 6 is dimmest to human naked eye.",          0, "15.1");
	options.Add("maxstarsize",  's', 1, "Pixels wide of brightest star",0, "10");
	options.Add("maxstarangle", 'a', 1, "Arc minutes wide of brightest star",0, "20");
	options.Add("bigthreshhold",'S', 1, "Magnitude brighter than this gets special size treatment",    0, "4");
	options.Add("alpha-amp",    'A', 1, "Adjustment to artifically pump up dim stars. Between 0 and 1",0, "0");
	options.Add("halo",         'H', 1, "Render halos on bright stars, number is ratio of "
		                                "halo diameter to star diameter",                              0, "7");
	options.Add("transparent",  'T', 0, "Put stars on transparency rather than a black background",    0, NULL);
	options.Add("galactic",     'g', 0, "Make it so the Milky Way is horizontal",                      0, NULL);

	options.Add("help",         'h', 0, "Show help and exit",                                          0, NULL);
	options.Add("help-html",    'l', 0, "Output help in html form and exit",                           0, NULL);
	options.Add("version",      'v', 0, "Show version and exit",                                       0, NULL);
}


int main(int argc, char **argv)
{
	InitializeMagick(*argv);


	RenderContext rr;

	 //process command line options
	InitOptions();
	int cc,index;

	cc=options.Parse(argc,argv, &index);
    if (cc==-2) {
        cerr <<"Missing parameter for "<<argv[index]<<"!!"<<endl;
        exit(0);
    }
    if (cc==-1) {
        cerr <<"Unknown option "<<argv[index]<<"!!"<<endl;
        exit(0);
    }


	 //initialize option variables
	const char *tycho_file=NULL;
	const char *pgc_file=NULL;
	double ang=-1;
	double alphaamp=0;


    LaxOption *o;
    for (o=options.start(); o; o=options.next()) {
        switch(o->chr()) {
            case 'l': {
				options.HelpHtml(stdout);
				exit(0);
			  }
            case 'h': {
				options.Help(stdout);
				exit(0);
			  }
            case 'v': {
				cout << options.HelpHeader()<<endl;
				exit(0);
			  }

            case 'w': {
				rr.width=strtol(o->arg(),NULL,10);
				if (rr.width<2) { cerr <<"Badth width value!"<<endl; exit(1); }
				rr.height=rr.width/2;
              } break;

            case 'o': {
				makestr(rr.filename,o->arg());
              } break;

            case 'p': {
				pgc_file=o->arg();
              } break;

            case 'c': {
				tycho_file=o->arg();
              } break;

            case 'm': {
				rr.minimum_magnitude=strtod(o->arg(),NULL);
              } break;

            case 'g': {
				rr.galactic=1;
              } break;

            case 'T': {
				rr.transparent=1;
              } break;

            case 'H': {
				rr.usehalo=strtod(o->arg(),NULL);
              } break;

            case 's': {
				rr.maxstarsize=strtod(o->arg(),NULL);
              } break;

            case 'a': {
				ang=strtod(o->arg(),NULL);
              } break;

            case 'S': {
				rr.bigthreshhold=strtod(o->arg(),NULL);
              } break;

            case 'A': {
				alphaamp=strtod(o->arg(),NULL);
				cerr <<" *** need to implement putting on base alphaamp to point amp curve:"<<alphaamp<<endl;
              } break;

		}
	}

	for (o=options.remaining(); o; o=options.next()) {
		FILE *f=fopen(o->arg(),"r");
		if (f) {
			rr.dump_in(f,0,0,NULL,NULL);
			fclose(f);
		}
	}


	if (ang>0) {
		rr.maxstarsize=rr.width*ang/60./360.;
	}
	if (rr.usehalo>1) {
		rr.halowidth=rr.maxstarsize*rr.usehalo;
	}


	//pixelwidth=1./width;



	 //-------print options summary
	cout <<endl;
	cout << "      Tycho catalog: "<<(tycho_file?tycho_file:"none")<<endl;
	cout << "     Galaxy catalog: "<<(pgc_file?pgc_file:"none")<<endl;
	cout << "      Outputting to: "<<rr.filename<<endl;
	cout << "         Dimensions: "<<rr.width<<" x "<<rr.height<<endl;
	cout << "Brightness at least: "<<rr.minimum_magnitude<<endl;
	cout << "        Halos after: "<<rr.bigthreshhold<<endl;


	 //fill in up render context
	if (pgc_file)   rr.catalogs.push(new Catalog("Principal Galaxy Catalog", pgc_file,   PrincipalGalaxy),1);
	if (tycho_file) rr.catalogs.push(new Catalog("Tycho 2 Star Catalog",     tycho_file, Tycho2),         1);

	//rr.catalogs.push(new Catalog("Custom Galaxy",     "galaxy.dat", CustomGalaxy), 1);

	 //------set up windows...
	anXApp app;
	makestr(app.app_profile,"Dark");
	app.init(argc,argv);

	app.addwindow(new MainWindow(&rr));
	

	 //...and off we go!
	app.run();


	cout <<"--- Shutting down... ----\n";
	app.close();
	
	cout <<"------ Bye! -------\n";
	return 0;

} //main()










