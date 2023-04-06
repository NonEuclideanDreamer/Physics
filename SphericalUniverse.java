
	import java.awt.Color;
	import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
	import java.util.Random;

import javax.imageio.ImageIO;
 
public class SphericalUniverse 
{
	public static int blue=Color.black.getRGB(),
					gray=Color.darkGray.getRGB(),  
			hundred=255; 
	public static boolean relativistic=false,
				manual=true;
	public static String name="spher";
	public static int background=new Color(0,0,0).getRGB(); 
	public static int[] bg= {0,0,0},  
			bl= {128,128,128}; 
	public static int n=5, //number of objects 
				dim=2,
				counter=1, //starting time
				T=878, //iterations 
				width=1440,height=1440, 
				ten=10000;//what mass is considered considerable 
	static double scale=2*Math.PI; 
	public static double t=1, //time step
				blur=0.9,//0:only current position visible, 1: whole trace visible
				rulescope=1, 
				c=0.01,//speed of light 
				g=0.0001;//gravitational constant 
  
	//starting condition for manually set particles 279
	public static double[][]loc0={{0,0.1},{-Math.PI/3,Math.PI/2-0.5},{Math.PI,Math.PI/2-0.5},{0,Math.PI*5/4},{2,0},{42,0},{63,0}},//,{0,50},{0,10}},//
				inertia0={{0,0},{0.1,0.02},{0.1,0.02},{0,-0.1},{-0,0.1},{0,-0.1},{0,0.1}} ;
	public static double[] mass= {1,1,1,1,1,1,1},
				radius= {0.01,0.01,0.01,0.01,0.01,1,1};
	public static int[][] color= {{128,128,128},{255,0,0},{0,255,0},{0,0,255},{255,255,255},{255,0,255},{255,255,0},{0,255,255}};
	public static BufferedImage screen=new BufferedImage(width,height,BufferedImage.TYPE_4BYTE_ABGR);

	public static void main(String[] args) 
	{ 
		//Particle.
		Cluster.fractal=false;
		Cluster.rulescope=rulescope;
		Cluster.relativistic=relativistic;
		Particle.spherical=true;
		Particle.r=1;
		Screen.height=height;
		Screen.width=width;
		Screen.blur=blur;
		Random rand=new Random();
		Particle.g=g;Particle.c=c;
		
		Screen.scale=scale;
		ArrayList<Screen> screens=new ArrayList<Screen>();
		screens.add(new Screen(Particle.spherzero()));
		Cluster cluster;
		if(manual)
		{	ArrayList<Particle> obj=new ArrayList<Particle>();
		for(int i=0;i<n;i++) 
		{

			obj.add(new Particle(loc0[i],inertia0[i],mass[i],radius[i],color[i]));
			
		} 
		 
		cluster=new Cluster(obj,screens);
		}
		else
		{
			double vmid=0.025,vsd=0.025,mass=0.00001; 
			double[]zro= {0,0},highloc= {scale/2,scale/2},highv= {1,1};
			Particle mid=new Particle(zro,zro,0.5,0.5,background),high=new Particle(highloc,highv,0.5,0.5,background) ,center=new Particle(zro,zro,1,1,background);
			cluster=Cluster.spherbigBang(n,vmid,vsd,mass,name);//Cluster.spherCloud(n);//	
			//random: Particle cloud vs bigBang: Start in the middle with outwards velocity
		}
		
		{
			
		cluster.screens.get(0).focus.name=name;
	
		int[][][]front	=new int[width][height][2],
				back=new int[width][height][2];
	cluster.print();
	double truetime=0.1,minrad=Math.PI;int steps;
		for(int time=counter;time<T+counter;time++)//for(int k=0;k<7;k++) 
		{
			//Particle.r=12*(Math.sin(time/100.0));
			//cluster.adjustRadii();
			//Particle.g=Math.pow((24.0/time),2);
			//Particle.c=24.0/time;
		//	System.out.println(time+",Universesize="+Particle.r+" n="+cluster.objects.size()+", x="+cluster.objects.get(0).loc[0]+", y="+cluster.objects.get(0).loc[1]);
			//
		
		if(time>438)			draw(cluster,front,back,time); 
			//drawForcefield(cluster, time); 
				//c/R*step<minrad->step<minrad*R/c ->t/step>t*c/minrad/R
			/*	steps=1+(int)(t*c/minrad/Particle.r);
				System.out.println("steps"+steps);
				for(int i=0;i<steps;i++)*/
		//cluster.printcomplexcoord();
			cluster.sphericalupdate(t,Particle.r); 
		//	minrad=cluster.minrad();
			
			
			
		}
		}
	} 
	
	public static void drawForcefield(Cluster cluster, int time)
	{
		int size=1440,
				n=cluster.objects.size();
		BufferedImage image=new BufferedImage(size,size,BufferedImage.TYPE_4BYTE_ABGR);
		double trans=0.8;
		
		
		Particle probe=Particle.zro(), probe2=Particle.zro();

		for(int i=0;i<size;i++)
		for(int j=0;j<size;j++)
		{	double rad=(Math.pow(i-size/2, 2)+Math.pow(j-size/2, 2))/size/size*4;
			if(rad<=1)
		{
			double[]front=new double[2],back=new double[2];//color/opaquicity
			/*probe.loc[0]=Math.atan2(j-size/2, i-size/2);
			probe.loc[1]=Math.asin(Math.sqrt(rad));
			probe2.loc[0]=probe.loc[0];
			probe2.loc[1]=Math.PI-probe.loc[1];*/
			probe.loc[1]=Math.acos((2.0*j-size)/size);
			probe.loc[0]=Math.acos((2.0*i-size)/size/(Math.sqrt(1-Math.pow(2.0*j/size-1, 2))));
			probe2.loc[0]=2*Math.PI-probe.loc[0];
			probe2.loc[1]=probe.loc[1];
			double[]direction=new double[dim];//for torus 
	
			for(int k=0;k<n;k++)
			{
				Particle object=cluster.objects.get(k);
			
				
				for(int l=0;l<dim;l++)
				{
					direction[l]=probe.loc[l]-object.loc[l];
				}
				
				double[]acc	= object.spherAcc(probe, direction, false, 1);
				for(int l=0;l<2;l++)
					front[l]+=acc[l]*object.mass;
				
				for(int l=0;l<dim;l++)
				{
					direction[l]=probe2.loc[l]-object.loc[l];
				}
				
				acc	= object.spherAcc(	probe2, direction, false, 1);
				for(int l=0;l<2;l++)
					back[l]+=acc[l]*object.mass;
			}
			image.setRGB(i,j,color(huecolor(front),huecolor(back),rad));
		}
		else image.setRGB(i, j, gray);
		}
		File outputfile = new File("spherForce"+time+".png");
		try 
		{
			ImageIO.write(image, "png", outputfile);
		} catch (IOException e) {
			System.out.println("IOException");
			e.printStackTrace(); 
		}	
	}
	
	public static void draw(Cluster cluster,int[][][]front, int[][][]back,int time)
	{
		int size=height,
				n=cluster.objects.size(),
				border=(width-size)/2,
				yb=(height-size)/2;
		double trans=0.8;
		
		for(int i=0;i<size;i++)
		{
			for(int j=0;j<size;j++)
			{
				front[i][j][1]*=blur;
				back[i][j][1]*=blur;
			}
		}
		
		BufferedImage image=new BufferedImage(width,height,BufferedImage.TYPE_4BYTE_ABGR);
		
		for(int i=0;i<n;i++)
		{
	
			double[] loc=cluster.objects.get(i).loc,
					center= {Math.cos(loc[0])*Math.sin(loc[1]),Math.sin(loc[0])*Math.sin(loc[1])};
			double a=Math.pow(Math.cos(loc[0])/Math.cos(loc[1]),2)+Math.pow(Math.sin(loc[0]), 2),
					b=2*Math.sin(loc[0])*Math.cos(loc[0])*Math.pow(Math.tan(loc[1]), 2),
					c=Math.pow(Math.sin(loc[0])/Math.cos(loc[1]), 2)+Math.pow(Math.cos(loc[0]),2),
					r=cluster.objects.get(i).radius;
			int color=cluster.objects.get(i).color;
			//System.out.println("center="+center[0]+","+center[1]+", r="+r+", phi="+loc[0]+",psi="+loc[1]+", a="+a);
			for(int y=(int) (-size*r/2);y<size*r/2;y++)
			{
				//System.out.print("y="+y);
				double det=Math.pow(b*y, 2)-4*a*(c*y*y-r*r*size*size/4);
				if(det>0)
				{
					det=Math.sqrt(det)/2/a;
					double s=-b*y/a/2;
					//System.out.println("s="+s+", det="+det);
					for(int x=(int)(size*center[0]/2+s-det);x<size*center[0]/2+s+det;x++)
					{
						//System.out.println("x,y="+x+","+y);
						if(loc[1]<Math.PI/2)
						{try 
							{front[x+size/2][y+(int)((center[1]+1)*size/2)]=new int[] {color,hundred};}
						catch(ArrayIndexOutOfBoundsException e) {}
						}	
						else
							try{back[x+size/2][(int) (y+(center[1]+1)*size/2)]=new int[] {color,hundred};}
					
					catch(ArrayIndexOutOfBoundsException e) {}
				}}
			}
		}
		
		for(int x=0;x<width;x++)
			if(x>(width-size)/2&&x<(width+size)/2)
			for(int y=0;y<height;y++)
			{
				double rad=(Math.pow(width/2-x, 2)+Math.pow(height/2-y, 2))/size/size*4;
				if(rad>1)image.setRGB(x,y,gray);
				else image.setRGB(x,y,color(front[x-border][y],back[x-border][y],rad));
			}
			else
				for(int y=0;y<height;y++)
				{
					image.setRGB(x,y,gray);
				}
		File outputfile = new File(name+time+".png");
		try 
		{
			ImageIO.write(image, "png", outputfile);
		} catch (IOException e) {
			System.out.println("IOException");
			e.printStackTrace(); 
		}	
	}

	private static int color(int[] is, int[] is2, double rad) 
	{
		Color col1=new Color(is[0]), 
				col2=new Color(is2[0]);
	
		int[] c1= {col1.getRed(),col1.getGreen(),col1.getBlue()},
				c2= {col2.getRed(),col2.getGreen(),col2.getBlue()},
				out=new int[3];
		double opaquicity=is[1]*(1-rad)+hundred*rad;
		for(int i=0;i<3;i++)
			out[i]=(int)(is[1]*c1[i]+(hundred-is[1])*(rad*bl[i]+(1-rad)*(is2[1]*c2[i]+(hundred-is2[1])*bg[i])/hundred))/hundred;
		//System.out.print("blue="+out[2]);
		return new Color(out[0],out[1],out[2]).getRGB();
	}
	public static int[] huecolor(double[] v)
	{
		double h=(Math.atan2(v[1], v[0])+Math.PI)/Math.PI*3,
				l=Math.sqrt(Math.pow(v[0],2)+Math.pow(v[1], 2))*100;
		int k=0;
		if (l>1)l=1;
		int[] c=new int[3];
		while(h>2)
		{
			h-=2;
			k++;
		}
		int m=k+1;
		if(m==3)m=0;
		if(h<1)
		{
			c[k]=(int)(255);
			c[m]=(int)(256*h);
		}
		else if(h==1)
		{
			c[k]=(int)(255);
			c[m]=(int)(255);
		}
		else
		{
			c[k]=(int)(256*(2-h));
			c[m]=(int)(255);
		}
		return new int[] {new Color(c[0],c[1],c[2]).getRGB(),(int)(l*hundred)};
	}
	
}


