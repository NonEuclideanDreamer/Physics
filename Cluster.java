import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import javax.imageio.ImageIO;

//**********************************************************************************************
//Author: Non-Euclidean Dreamer
//A Cluster of Particles acting on each other
//**********************************************************************************************


public class Cluster 
{
	ArrayList<Particle> objects;
	ArrayList<Screen> screens;
	int dim;
	static double rulescope=1;
	static double metric=1;
	static int black=Color.black.getRGB(), matters=0;	
	static boolean relativistic=true;
	public Cluster(ArrayList<Particle> obj)
	{
		objects=obj;
		dim=obj.get(0).dim;
		screens=new ArrayList<Screen>();
		screens.add(new Screen(Particle.zro()));
	}
	
	public Cluster(ArrayList<Particle> obj,ArrayList<Screen>scr)
	{
		objects=obj;
		dim=obj.get(0).dim;
		screens=scr;
	}
	
	public static Cluster regular(int hor,int ver)
	{
		ArrayList<Particle>obj=new ArrayList<Particle>();
		for(int i=0;i<hor;i++)
			for(int j=0;j<ver;j++)
			{
				double[] loc= {0+i*1000/hor,0+j*1000/ver},v= {0,0},loc1= {100+i*1000/hor,0+j*1000/ver},v1= {0,i%2*2-1};
				obj.add(new Particle(loc,v,((i+j)%2-0.5)*4,10,new Color((i+j)%2*255,255-(i+j)%2*255,0).getRGB()));
				obj.add(new Particle(loc1,v1,0,2,new Color(0,0,255).getRGB()));
			}
		return new Cluster(obj);
	}


	

	public static Cluster random(int n,int dim,Particle mid,Particle high)
	{
		Random rand=new Random();
		
		ArrayList<Particle> obj=new ArrayList<Particle>();
		for(int i=0;i<n;i++)
		{
			double[]loc=new double[dim],v=new double[dim];
			for(int j=0;j<dim;j++)
			{
				loc[j]=randomDouble(rand,mid.loc[j],high.loc[j]);
				v[j]=randomDouble(rand,mid.v[j],high.v[j]);
			}
			obj.add(new Particle(loc,v,Math.signum(randomDouble(rand,mid.mass,high.mass)),randomDouble(rand,mid.radius,high.radius),randomColor(rand)));
		}
		return new Cluster(obj);
	}
	
	public void torusupdate(double t,double[]size)
	{	
		int n=objects.size();
		Particle[] out=new Particle[n];
		
		for(int i=0;i<n;i++)
		{
			out[i]=objects.get(i).copy();
	
		}
		for(int i=0;i<n;i++)
		{
			double[] acc=new double[dim];
			for(int j=0;j<n;j++)
				if(i!=j&&out[j].mass!=0)
				{
					double d;
					double[]direction=new double[dim];//for torus
					for(int k=0;k<dim;k++)
					{
						direction[k]=out[i].loc[k]-out[j].loc[k];
						if(direction[k]>size[k]/2)
							direction[k]-=size[k];
						else if(direction[k]<=-size[k]/2)
							direction[k]+=size[k];
					}
					if(matters==0)d=out[i].torusdistance(direction, metric,size);
					else d=out[i].vdifference(out[j], metric);
					add(acc,out[i].torusgravAcc(out[j],d,direction,relativistic));
				
				}
			if(relativistic)objects.get(i).relPull(acc, t);
			else objects.get(i).pull(acc,t);
			//objects.get(i).shove(acc,t);
		}    
		
		for(int i=0;i<n;i++)
		{
			boolean delete=false;
			for(int j=objects.size()-1;j>i;j--)
				if(objects.get(i).distance(objects.get(j), metric)<objects.get(i).radius+objects.get(j).radius) 
				{ 
					delete=objects.get(i).merge(objects.get(j),this,delete);
					objects.remove(j);	
				}
			if(delete) {objects.remove(i);i--;}
		}
		sort();
		saturize();

		System.out.print("biggest="+objects.get(objects.size()-1).name+" with mass="+objects.get(objects.size()-1).mass);
	}
	public int update(double t)
	{	
		int n=objects.size();
		Particle[] out=new Particle[n];
		
		for(int i=0;i<n;i++)
		{
			out[i]=objects.get(i).copy();
	
		}
		for(int i=0;i<n;i++)
		{
			if(relativistic)
			{
				double[]force=new double[dim];
				for(int j=0;j<n;j++)
					if(i!=j&&out[j].mass!=0)
					{
						add(force,out[i].relgravForce(out[j]));
					}
				objects.get(i).relPull(force,t);
			}
			else
			{
			double[] acc=new double[dim];
			for(int j=0;j<n;j++)
				if(i!=j&&out[j].mass!=0)
				{
					add(acc,out[i].gravAcc(out[j]));
				}
			objects.get(i).pull(acc,t);
			//objects.get(i).shove(acc,t);
			}
		}    
		
		for(int i=0;i<n;i++)
		{
			boolean delete=false;
			for(int j=objects.size()-1;j>i;j--)
				if(objects.get(i).distance(objects.get(j), metric)<objects.get(i).radius+objects.get(j).radius) 
				{ 
					if(i==0)return objects.get(j).color;//For Fractal
					delete=objects.get(i).merge(objects.get(j),this,delete);
					objects.remove(j);	
				}
			if(delete) {objects.remove(i);i--;}
		}
        sort();
		saturize();

		//System.out.print("biggest="+objects.get(objects.size()-1).name+" with mass="+objects.get(objects.size()-1).mass);
		return 0;
	}


	private void add(double[] d1, double[] d2) 
	{
		for (int i=0;i<d1.length;i++)
		{
			d1[i]+=d2[i];
		}
		
	}

	private void saturize() 
	{
		for(int i=0;i<objects.size();i++)
		{
			objects.get(i).saturize();
		}
		
	}

	public void draw(BufferedImage screen, int scale, int background, String name,int t) 
	{
		int height=screen.getHeight(),width=screen.getWidth();
		for(int i=0;i<width;i++)
		{
			for(int j=0;j<height;j++)
			{
				screen.setRGB(i, j, background);
			}
		}
		for(int k=0;k<objects.size();k++)
		{	Particle obj= objects.get(k);
			double[]center= {width/2+obj.loc[0]/scale*width,height/2+obj.loc[1]/scale*width};
			double r=obj.radius*width/scale;
			for(int i=(int)(center[0]-r);i<center[0]+r+1;i++)
			{
				double s=Math.pow(Math.pow(r, metric)-Math.pow(Math.abs(center[0]-i), metric), 1.0/metric);
				for(int j=(int)(center[1]-s);j<center[1]+s+1;j++)
				{try {
					screen.setRGB(i,j,obj.color);}catch(ArrayIndexOutOfBoundsException e) {}
				}
			}
				
		}                                                   
		File outputfile = new File(name+t+".png");
		try 
		{
			ImageIO.write(screen, "png", outputfile);
		} catch (IOException e) {
			System.out.println("IOException");
			e.printStackTrace(); 
		}
	}


	public void draw( int t) 
	{
		for(int i=0;i<screens.size();i++)
		screens.get(i).draw(this,t);
		
	}


	public static int blur(int col1, int col2,double blur) 
	{
		Color c1=new Color( col1), c2=new Color(col2);
		int red=(int)(c1.getRed()*blur+(1-blur)*c2.getRed()),blue=(int)(c1.getBlue()*blur+(1-blur)*c2.getBlue()),green=(int)(c1.getGreen()*blur+(1-blur)*c2.getGreen());
		return new Color(red,green,blue).getRGB();
		
	}


	public void centralize() 
	{
		for(int i=0;i<dim;i++)
		{
			double x=0,v=0,m=0;
			for(int j=0;j<objects.size();j++)
			{
				Particle obj=objects.get(j);	
				double mass=Math.abs(obj.mass);
				x+=mass*obj.loc[i];
				v+=mass*obj.v[i];
				m+=mass;
			}
			x/=m;v/=m;System.out.println("shift="+x+", move="+v+", total mass="+m);
			for(int j=0;j<objects.size();j++)
			{
				Particle obj=objects.get(j);
				obj.loc[i]-=x;
				obj.v[i]-=v;
			}
		}
		
	}
	public static double randomDouble(Random rand,double mid,double high)
	{
		return (rand.nextDouble()-0.5)*2*(high-mid)+mid;
	}

	//Gives saturated color
	public static int randomColor(Random rand)
	{
		double red=rand.nextDouble()-0.5,green=rand.nextDouble()-0.5, blue=rand.nextDouble()-0.5,max=Math.max(Math.abs(red), Math.max(Math.abs(green), Math.abs(blue)));
		 red=red*127.9/max+128;green=green*127.9/max+128;blue=blue*127.9/max+128;
		 return new Color((int)red,(int)green,(int)blue).getRGB();
	}

	public static Cluster bigBang(int n, double vmid, double vsd, String name) 
	{
		ArrayList<Particle> obj=new ArrayList<Particle>();
		Random rand=new Random();double[]loc= {0,0};
		for (int i=0;i<n;i++)
		{		
			double angle=rand.nextDouble()*2*Math.PI,norm=rand.nextGaussian(vmid, vsd);
			double[]v={Math.cos(angle)*norm,Math.sin(angle)*norm},lo= {2*v[0],2*v[1]};
			obj.add(new Particle(lo,v,1,1,new Color(255*rand.nextInt(2),255*rand.nextInt(2),255*rand.nextInt(2)).getRGB()));//new Color(128+(int)(128*Math.cos(angle)),128+(int)(128*Math.sin(angle)),Math.max(Math.min(255,128+(int)((norm-vmid)*128/3/vsd)),0)).getRGB()));
		}
		return new Cluster(obj);
	}
	public void sort()
	{
		Collections.sort(objects);
	}

	public void replaceScreen(Particle object,Particle nju) 
	{
		boolean notfound=true;
		int i=0;
		
		while(notfound)
		{
			
			if(screens.get(i).focus==object) 
			{
				screens.add(new Screen(screens.get(i).scr,nju));
				screens.remove(i);
				return;
			}
			else i++;
		}
		
	}
	public void print()
	{
		for(int i=0;i<objects.size();i++)
			objects.get(i).print();
	}
	public void torus(double[]size)
	{
		double m=size[0],n=size[1];
		double[] loc;
		for(int i=0;i<objects.size();i++)
		{
			loc=objects.get(i).loc;
			loc[0]=((loc[0]+m/2)%m+m)%m-m/2;
			loc[1]=((loc[1]+n/2)%n+n)%n-n/2;
		
		}
	}
	
	public void updateCharged(double t)
	{
	int n=objects.size();
	Particle[] out=new Particle[n];
	
	for(int i=0;i<n;i++)
	{
		out[i]=objects.get(i).copy();

	}
	for(int i=0;i<n;i++)
	{
		double[] acc=new double[dim];
		for(int j=0;j<n;j++)
			if(i!=j&&out[j].mass!=0)
			{
				add(acc,out[i].gravAcc(out[j]));
				add(acc,(out[i]).elAcc(out[j]));
			}
		objects.get(i).pull(acc,t);
		//objects.get(i).shove(acc,t);
	}    
	
	for(int i=0;i<n;i++)
	{
		boolean delete=false;
		for(int j=objects.size()-1;j>i;j--)
			if(objects.get(i).distance(objects.get(j), metric)<objects.get(i).radius+objects.get(j).radius) 
			{ 
				delete=objects.get(i).merge(objects.get(j),this,delete);
				objects.remove(j);	
			}
		if(delete) {objects.remove(i);i--;}
	}
	sort();
	saturize();

	System.out.print("biggest="+objects.get(objects.size()-1).name+" with mass="+objects.get(objects.size()-1).mass);
}

	public static Cluster randomch(int n, int dim, Particle mid, Particle high) 
	{

		Random rand=new Random();
		
		ArrayList<Particle> obj=new ArrayList<Particle>();
		for(int i=0;i<n;i++)
		{
			double[]loc=new double[dim],v=new double[dim];
			for(int j=0;j<dim;j++)
			{
				loc[j]=randomDouble(rand,0,high.loc[j]);
				v[j]=randomDouble(rand,0,high.v[j]);
			}Particle o=new Particle(loc,v,randomDouble(rand,mid.mass,high.mass),randomDouble(rand,mid.radius,high.radius),randomColor(rand));
			o.setCharge(rand.nextInt(3)-1);
			obj.add(o);
		}
		return new Cluster(obj);
	}

	public boolean escapes(int i)
	{
		boolean out=false;
		return out;
	}
}
