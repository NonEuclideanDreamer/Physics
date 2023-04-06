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
	public double mass=5;//Gesamtmasse
	public static boolean[]mirror= {false,false};//For toroidal{0,0},for Klein {1,0} for RP2 {1,1}
	static double rulescope=1;
	static double metric=2;
	static int black=Color.black.getRGB(), matters=0;	
	public static boolean relativistic=true,
			fractal;
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
	public static Cluster spherCloud(int n)
	{
		Random rand=new Random();
		
		ArrayList<Particle> obj=new ArrayList<Particle>();
		double m=0.0001, r=r(0,m,1);
		for(int i=0;i<n;i++)
		{
			double theta=Math.acos(1-2*rand.nextDouble());
			double phi=(rand.nextDouble()*2-1)*Math.PI;
			
			obj.add(new Particle(new double[] {phi,theta},new double[] {0.01,0},m,r,new int[] {rand.nextInt(2)*255,rand.nextInt(2)*255,rand.nextInt(2)*255}));
		}
		return new Cluster(obj);
	}
	public void update(double t, Gradient gr)
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
					double[]r=out[j].direction(out[i]);
					
					double[][] rep=representant(gr,r);
					r=rep[0];
					double[]sign=rep[1];//are we turning the gradient around?
					//System.out.println("r=("+r[0]+","+r[1]+")");
					//System.out.println("sign=("+sign[0]+","+sign[1]+")");
					double[]d=gr.grad[(int)Math.round(r[0]*gr.res)][(int)Math.round(r[1]*gr.res)];
					
					for(int k=0;k<dim;k++)
					{
						acc[k]+=d[k]*sign[k]*out[j].mass*Particle.g;
					}
						//System.out.println("What?");
				
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

	}
	private double[][] representant(Gradient gr,double[] r) 
	{
		double[][]out={{r[0],r[1]}, {1,1}};
	
	
		for(int k=0;k<dim;k++)
		{
			if(out[0][k]<0) {out[0][k]*=-1; out[1][k]*=-1;}
		
			
		}
		return out;
	}

	private static double[] times(double[] ds, double m) 
	{
		double[]out=new double[ds.length];
		for(int i=0;i<ds.length;i++)
			out[i]=ds[i]*m;
		return out;
	}
	public int sphericalupdate(double t,double r)
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
			acc[0]=0;acc[1]=0;
			for(int j=0;j<n;j++)
				if(i!=j&&out[j].mass!=0)
				{
					double[]direction=new double[dim];//for torus
					for(int k=0;k<dim;k++)
					{
						direction[k]=out[j].loc[k]-out[i].loc[k];
					}
					
					add(acc,out[i].spherAcc(out[j],direction,relativistic,r));
				
				}
			if(relativistic) {objects.get(i).relPull(acc, t);System.out.println("relativistic sphere2 not yet implemented");}
			else objects.get(i).spherpull(acc,t,r);
			//objects.get(i).shove(acc,t);
		}    
		
		for(int i=0;i<n;i++)
		{
			boolean delete=false;
			for(int j=objects.size()-1;j>i;j--)
				if(objects.get(i).sphericalDistance(objects.get(j))<objects.get(i).radius+objects.get(j).radius) 
				{ 
					if(i==0&&fractal)return objects.get(j).color;//For Fractal
					delete=objects.get(i).merge(objects.get(j),this,delete);
						
					
					objects.remove(j);
				}
			if(delete) {objects.remove(i);i--;}
		}
		n=objects.size();
		Random rand=new Random();
		double c=Particle.c,R=Particle.r;
		for(int i=n-1;i>-1;i--)
		{
			double u=objects.get(i).u();
			//System.out.println("u="+u+", mR^2="+(mass*R*R));
			double random=rand.nextDouble();
			double stability=5.0/t;
			if(random>Math.exp(-u/stability*0))
			{	
				
				random=-Math.log(random)*stability;
				System.out.println("random="+random);
				//if(random>0)
				{
					
					double mass=objects.get(i).mass;
					double massratio=rand.nextDouble();// ToDo: Maybe needs to change
				double excess=u-random;//This is the new kinetic energy
				random=u-excess;
				System.out.println("new u="+random);	
				double m1=mass*massratio/(massratio+1),m2=mass-m1;
				double v1=Math.sqrt(excess*(2*c*c*m2+excess)*(2*c*c*m1+excess)*(2*c*c*mass+excess))/(2*c*m1*(c*c*mass+excess));
				System.out.println("v1="+v1);
				double angle=rand.nextDouble(2*Math.PI);
				if(Particle.norm(objects.get(i).v)!=0) {double sign=Math.signum(angle-Math.PI);angle=Math.atan2(sign*objects.get(i).v[1],-sign*objects.get(i).v[0]);}
				double[]vel1= {Math.cos(angle)*v1,Math.sin(angle)*v1},
						vel2= {-massratio*vel1[0],-massratio*vel1[1]};
				double u2=random/(1+massratio),
						u1=u2*massratio,
						r1=r(u1,m1,objects.get(i).radius),
						r2=r(u2,m2,objects.get(i).radius);
				System.out.println("u1="+u1+", u2="+u2+",r1="+r1+", r2="+r2);
				double[]vl1=objects.get(i).v.clone(),vl2=objects.get(i).v, loc=objects.get(i).loc;
				add(vl1,vel1);add(vl2,vel2);
				int[]co=colorsplit(objects.get(i).color,massratio/(massratio+1));
				Particle p1 =new Particle(loc,vl1,m1,r1,co[0]);
				Particle p2 =new Particle(loc,vl2,m2,r2,co[1]);
				int d=0;
				while(p1.sphericalDistance(p2)<r1+r2&&d<1000)
				{
					p1.spherpush(vel1);
					p2.spherpush(vel2);
					d++;
				}
				
				objects.add(p1);
				objects.add(p2);
				
				objects.remove(i);
			}
		}}
		
		if(!fractal)sort();
		saturize();
	
		//System.out.print("biggest="+objects.get(objects.size()-1).name+" with mass="+objects.get(objects.size()-1).mass);
		return 0;
	}
	private int[] colorsplit(int color, double mass1) 
	{
		Random rand=new Random();
		double down=rand.nextDouble();
		int[][]out=new int[2][3];
		Color c=new Color(color);
		double[]cl= { (c.getRed()-128)*down,
			(c.getGreen()-128)*down,
			(c.getBlue()-128)*down};
		double[]mass= {mass1,1-mass1};
		
		
		
		for(int i=0;i<3;i++)
		{
			int k=rand.nextInt(2);
			boolean b=rand.nextBoolean();
			
			if(b) 
			{
				out[k][i]=Math.min(127, (int)((cl[i]+mass[1-k]*128)/mass[k]));
			}
			else
			{
				out[k][i]=Math.max(-128, (int)((cl[i]-mass[1-k]*127)/mass[k]));
			}				
			out[1-k][i]=(int) ((cl[i]-mass[k]*out[k][i])/mass[1-k]);
		}
		int[]o=new int[2];
		for (int i=0;i<2;i++)
		{
			//System.out.print("blue="+out[i][2]);
			Color clr=new Color(Math.min(255,Math.max(0,out[i][0]+128)),Math.min(255,Math.max(0,out[i][1]+128)),Math.min(255,Math.max(0,out[i][2]+128)));
			o[i]=clr.getRGB();
		}
		return o;
	}

	public static double r(double u1, double m1,double limit)
	{
		double r=0,R=Particle.r,step=0.01,u;
		while(step>0.0000000000001)
		{
			
			do {
					r+=step;
					u=m1*(R*R+m1*Particle.g/2*(1+Math.cos(r))*(Math.log(1-Math.cos(r))-(1+Math.cos(r))/(1-Math.cos(r))));
					}
			while(u<u1&&r<limit);

			r-=step;
			step*=0.1;
		}

		return r;
	}
	public int finiteUpdate(double t, Potential p, Potential change)
	{
		int n=objects.size();
		Particle[] out=new Particle[n];
		
		for(int i=0;i<n;i++)
		{
			out[i]=objects.get(i).copy();
	
		}
		for(int i=0;i<n;i++)
		{
			
			double[]force=new double[2];
			int x=(int)Math.round(p.res*out[i].loc[0]),
					y=(int)Math.round(p.res*out[i].loc[1]);
			force[0]=-(p.right(x,y)-p.left(x,y))/2*p.res*mass;
			force[1]=-(p.up(x,y)-p.down(x,y))/2*p.res*mass;
			if(relativistic)objects.get(i).relPull(force, t);
			else objects.get(i).pull(force,t);
			//objects.get(i).shove(acc,t);
			objects.get(i).addPotential(change);
			System.out.println(p.right(x, y)+","+p.left(x, y));
		}    
		
		for(int i=0;i<n;i++)
		{
			boolean delete=false;
			for(int j=objects.size()-1;j>i;j--)
				if(objects.get(i).distance(objects.get(j), metric)<objects.get(i).radius+objects.get(j).radius) 
				{ 
					if(i==0&&fractal)return objects.get(j).color;//For Fractal
					delete=objects.get(i).merge(objects.get(j),this,delete);
						
					
					objects.remove(j);
				}
			if(delete) {objects.remove(i);i--;}
		}
		
	/*	n=objects.size();
		Random rand=new Random();
		for(int i=0;i<n;i++)
		{
			double u=objects.get(i).u();
			double random=rand.nextDouble();
			if(random<Math.exp(u))
			{
				double excess=u-Math.log(random);//This is the new kinetic energy
			}
		}*/
		if(!fractal)sort();
		saturize();
	
		//System.out.print("biggest="+objects.get(objects.size()-1).name+" with mass="+objects.get(objects.size()-1).mass);
		return 0;
	}
	public int torusupdate(double t,int[] torussize)
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
					double[]direction=new double[dim];//for torus
					for(int k=0;k<dim;k++)
					{
						direction[k]=out[j].loc[k]-out[i].loc[k];
					}
					
					add(acc,Particle.coverAcc(out[j],direction,relativistic,torussize,mirror,out[i].mass));
				
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
					if(i==0&&fractal)return objects.get(j).color;//For Fractal
					delete=objects.get(i).merge(objects.get(j),this,delete);
						
					
					objects.remove(j);
				}
			if(delete) {objects.remove(i);i--;}
		}
	/*	n=objects.size();
		Random rand=new Random();
		for(int i=0;i<n;i++)
		{
			double u=objects.get(i).u();
			double random=rand.nextDouble();
			if(random<Math.exp(u))
			{
				double excess=u-Math.log(random);//This is the new kinetic energy
			}
		}*/
		if(!fractal)sort();
		saturize();
	
		//System.out.print("biggest="+objects.get(objects.size()-1).name+" with mass="+objects.get(objects.size()-1).mass);
		return 0;
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
					if(i==0&&fractal)return objects.get(j).color;//For Fractal
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
		int red=(int)Math.round(c1.getRed()*blur+(1-blur)*c2.getRed()),blue=(int)Math.round(c1.getBlue()*blur+(1-blur)*c2.getBlue()),green=(int)Math.round(c1.getGreen()*blur+(1-blur)*c2.getGreen());
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
		int red=255*rand.nextInt(2),green=rand.nextInt(2)*255, blue=rand.nextInt(2)*255;
		 return new Color(red,green,blue).getRGB();
	}

	public static Cluster bigBang(int n, double vmid, double vsd, String name) 
	{
		ArrayList<Particle> obj=new ArrayList<Particle>();
		Random rand=new Random();double[]loc= {0,0};
		for (int i=0;i<n;i++)
		{		
			double angle=rand.nextDouble()*2*Math.PI,norm=rand.nextGaussian(vmid, vsd);
			double[]v={Math.cos(angle)*norm,Math.sin(angle)*norm},lo= {10*v[0],10*v[1]};
			obj.add(new Particle(lo,v,0.1,1,new Color(255*rand.nextInt(2),255*rand.nextInt(2),255*rand.nextInt(2)).getRGB()));//new Color(128+(int)(128*Math.cos(angle)),128+(int)(128*Math.sin(angle)),Math.max(Math.min(255,128+(int)((norm-vmid)*128/3/vsd)),0)).getRGB()));
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
	public void printcomplexcoord()//Print the coefficients of the complex function belonging to this state.
	{
		int i=0;
		double[]z=new double[2];
		System.out.print("{");
		for(Particle obj: objects)
		{
			for(int j=0;j<obj.mass;j++)
			{
				z=obj.z();
				System.out.print("{"+z[0]+","+z[1]+"}, ");
				i++;
			}
			
		}
		System.out.println("},");
	}
	public void print()
	{
		for(int i=0;i<objects.size();i++)
			objects.get(i).print();
	}
	//puts the objects back into their screen depending on mirror
	public void torus(int[] torussize)
	{
		double[] loc;
		for(int i=0;i<objects.size();i++)
		{
			loc=objects.get(i).loc;
			int[]flip=new int[dim];
			for(int j=0;j<dim;j++)
			{
				if(mirror[j])
				{
					flip[1-j]=(int)Math.round(loc[j]/torussize[j])%2;
					
				}
				loc[j]=((loc[j]+torussize[j]/2)%torussize[j]+torussize[j])%torussize[j]-torussize[j]/2;

			}
			for(int j=0;j<dim;j++)
			if(flip[j]!=0)
			{
				loc[j]*=-1;
				objects.get(i).v[j]*=-1;
			}
		
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
	public void drawPotential(Potential p, int t)
	{
		Gradient g=new Gradient(p);
		BufferedImage scr=new BufferedImage(p.width*p.res,p.height*p.res, BufferedImage.TYPE_INT_RGB);
		double[][]field=new double[scr.getWidth()][scr.getHeight()];
		int n=objects.size();
		for(int k=0;k<n;k++)
		{
			double potential;
			Particle m=objects.get(k);
		for(int i=0;i<scr.getWidth();i++)
			for(int j=0;j<scr.getHeight();j++)
			{
				
				int[]index= {i,j};
				double[]direction=new double[2];
				for(int l=0;l<2;l++)
					direction[l]=m.loc[l]-index[l]*1.0/p.res+p.height/2;
				double[][]r=representant(g, direction);
				double d=p.p[(int)Math.round(r[0][0]*p.res)][(int)Math.round(r[0][1]*p.res)];
				
					field[i][j]+=d*m.mass;
			}
		}
		for(int i=0;i<scr.getWidth();i++)
			for(int j=0;j<scr.getHeight();j++)
			{
				
				int c=0;
				//	System.out.println(field[i][j]);
					//System.out.println(field[i][j][k]);
					if(field[i][j]<100)c=0;
					else if (field[i][j]>120)c=255;
					else c=(int)((field[i][j]-100)*256/20);		
			
				scr.setRGB(i,j,new Color(0, c, 0).getRGB());
			}
		File outputfile = new File("potentialfield"+t+".png");
		try 
		{
			ImageIO.write(scr, "png", outputfile);
		} catch (IOException e) {
			System.out.println("IOException");
			e.printStackTrace(); 
		}	
	}
	public void drawForce(Gradient p, int t)
	{
		BufferedImage scr=new BufferedImage(p.size[0]*p.res,p.size[1]*p.res, BufferedImage.TYPE_INT_RGB);
		double[][][]field=new double[scr.getWidth()][scr.getHeight()][2];
		int n=objects.size();
		for(int k=0;k<n;k++)
		{
			
			Particle m=objects.get(k);
		for(int i=0;i<scr.getWidth();i++)
			for(int j=0;j<scr.getHeight();j++)
			{
				
				int[]index= {i,j};
				double[]direction=new double[2];
				for(int l=0;l<2;l++)
					direction[l]=index[l]*1.0/p.res-p.size[l]/2-m.loc[l];
			//	double[][]r=representant(p,direction);
				double[]d=Particle.coverAcc(m, direction, relativistic, p.size, mirror,1);
				for(int l=0;l<2;l++)
					field[i][j][l]+=d[l]*m.mass;
			}
		}
		for(int i=0;i<scr.getWidth();i++)
			for(int j=0;j<scr.getHeight();j++)
			{
				double h=(Math.atan2(field[i][j][1], field[i][j][0])+Math.PI)/Math.PI*3,
						l=Math.sqrt(Math.pow(field[i][j][0],2)+Math.pow(field[i][j][1], 2));
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
					c[k]=(int)(255*l);
					c[m]=(int)(256*h*l);
				}
				else if(h==1)
				{
					c[k]=(int)(l*255);
					c[m]=(int)(l*255);
				}
				else
				{
					c[k]=(int)(256*l*(2-h));
					c[m]=(int)(l*255);
				}
				scr.setRGB(i,j,new Color(c[0],c[1],c[2]).getRGB());
			}
		File outputfile = new File("forcefield"+t+".png");
		try 
		{
			ImageIO.write(scr, "png", outputfile);
		} catch (IOException e) {
			System.out.println("IOException");
			e.printStackTrace(); 
		}	
	}

	public void drawForce(int t,double scale, int[]torussize, int mult) 
	{
		double br=5;//Brightness
		double[][][]field=new double[(int) (torussize[0]*scale)][(int) (torussize[1]*scale)][2];
		BufferedImage scr=new BufferedImage((int)(torussize[0]*scale*mult),(int)(torussize[1]*scale), BufferedImage.TYPE_INT_RGB);
		int n=objects.size();
		System.out.println("check1");
		for(int k=0;k<n;k++)
		{
			
			Particle m=objects.get(k);
		for(int i=0;i<scr.getWidth()/mult;i++)
			for(int j=0;j<scr.getHeight();j++)
			{
				
				int[]index= {i,j};
				double[]direction=new double[2];
				for(int l=0;l<2;l++)
					direction[l]=index[l]*1.0/scale-torussize[l]/2-m.loc[l];
			//	double[][]r=representant(p,direction);
				double[]d=Particle.coverAcc(m, direction, relativistic, torussize, mirror,1);
				for(int l=0;l<2;l++)
					field[i][j][l]+=d[l]*m.mass;
			}
		}
		for(int i=0;i<scr.getWidth()/mult;i++)
			for(int j=0;j<scr.getHeight();j++)
			{
				double h=(Math.atan2(field[i][j][1], field[i][j][0])+Math.PI)/Math.PI*3,
						l=br*Math.sqrt(Math.pow(field[i][j][0],2)+Math.pow(field[i][j][1], 2));
				int k=0;
				if (l>1)l=1;
				int[] c=new int[3],c2=new int[3];
				while(h>2)
				{
					h-=2;
					k++;
				}
				int m=k+1;
				if(m==3)m=0;
				if(h<1)
				{
					c[k]=(int)(255*l);
					c[m]=(int)(256*h*l);
				}
				else if(h==1)
				{
					c[k]=(int)(l*255);
					c[m]=(int)(l*255);
				}
				else
				{
					c[k]=(int)(256*l*(2-h));
					c[m]=(int)(l*255);
				}/*
				h=(Math.atan2(-field[i][j][1], field[i][j][0])+Math.PI)/Math.PI*3;
						l=br*Math.sqrt(Math.pow(field[i][j][0],2)+Math.pow(field[i][j][1], 2));
				k=0;
				if (l>1)l=1;
			
				while(h>2)
				{
					h-=2;
					k++;
				}
				 m=k+1;
				if(m==3)m=0;
				if(h<1)
				{
					c2[k]=(int)(255*l);
					c2[m]=(int)(256*h*l);
				}
				else if(h==1)
				{
					c2[k]=(int)(l*255);
					c2[m]=(int)(l*255);
				}
				else
				{
					c2[k]=(int)(256*l*(2-h));
					c2[m]=(int)(l*255);
				}*/
				scr.setRGB(i,j,new Color(c[0],c[1],c[2]).getRGB());
				//scr.setRGB(i+1*(int)(torussize[0]*scale),scr.getHeight()-1-j,new Color(c2[0],c2[1],c2[2]).getRGB());
			//	scr.setRGB(i+2*(int)(torussize[0]*scale), j, new Color(c[0],c[1],c[2]).getRGB());
			}
		File outputfile = new File("forcefield"+t+".png");
		try 
		{
			ImageIO.write(scr, "png", outputfile);
		} catch (IOException e) {
			System.out.println("IOException");
			e.printStackTrace(); 
		}	
	}

	public static Cluster spherbigBang(int n, double vmid, double vsd,double mass, String name) 
	{
		ArrayList<Particle> obj=new ArrayList<Particle>();
		Random rand=new Random();double[]loc= {0,0};
		double vangle=0.025,  r=r(0,mass,1);
		System.out.println("g="+Particle.g);
		for (int i=0;i<n;i++) 
		{		
			double angle=rand.nextDouble()*2*Math.PI,norm=rand.nextGaussian(vmid, vsd);
			double[]v={vangle,norm},lo= {angle,10*norm};
			obj.add(new Particle(lo,v,mass,r,new Color(255*rand.nextInt(2),255*rand.nextInt(2),255*rand.nextInt(2)).getRGB()));//new Color(128+(int)(128*Math.cos(angle)),128+(int)(128*Math.sin(angle)),Math.max(Math.min(255,128+(int)((norm-vmid)*128/3/vsd)),0)).getRGB()));
		}
		return new Cluster(obj);
	}
	public double minrad()//or avrad?
	{
		double out=Math.PI;
		for(Particle mass: objects)
		{
			out+=mass.radius;
		}
		return out/objects.size();
	}
	public void inflate(double r)
	{
		
	}

	public void adjustRadii() //If the universe shrinks, we have to make sure, the potential energy doesn't turn negative
	{
		for(Particle obj: objects)
		{
			if(obj.u()<0)
			{
				obj.radius=r(0,obj.mass,Math.PI);
			}
		}
		
	}
}
