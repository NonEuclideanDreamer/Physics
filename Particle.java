import java.awt.Color;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Random;

//**********************************************************************************************
// Author: Non-Euclidean Dreamer
// A physical Particle with all its basic properties
//**********************************************************************************************
public class Particle implements Comparable 
{
	static double g=6*Math.pow(10,-11),kc=9*Math.pow(10, 9),c=3*Math.pow(10, 4);
	double[]loc, v;//for relativistic universes v is used for p
	double mass, radius;
	static double rulescope=Cluster.rulescope, metric=1;
	int  color;
	static int dim=2,background=Color.black.getRGB();
	String name;
	Random rand=new Random();
	double charge;
	public Particle(double[]l,double[]vel,double m, double r,int c)
	{
		loc=l.clone();
		v=vel.clone();
		mass=m;
		radius=r;
		dim=l.length;
		color=c;
		name="";
	}
	public Particle copy()
	{
		Particle obj= new Particle(loc.clone(),v.clone(),mass, radius,color);
		obj.name=name.concat("");
		return obj;
	}
	public static Particle zro()
	{
		return new Particle(new double[dim],new double[dim],10,10,background);
	}
	public static Particle proton(double[]loc, double[]v)
	{
		Particle out=new Particle(loc,v,1,1,Color.red.getRGB());
		out.setCharge(1);
		return out;
	}
	public void setV(double[]vel)
	{
		v=vel.clone();
	}
	public double distance(Particle obj,double norm)
	{
		double out=0;
		for(int i=0;i<dim;i++)
		{
			out+=Math.pow(Math.abs(loc[i]-obj.loc[i]), norm);
		}
		out=Math.pow(out, 1/norm);
		
		return out;
	}
	public double vdifference(Particle obj, double norm)
	{
		double out=0;
		for(int i=0;i<dim;i++)
		{
			out+=Math.pow(Math.abs(v[i]-obj.v[i]), norm);
		}
		out=Math.pow(out, 1/norm);
		
		return out;
	}
	
	//Changes the inertia&location
	public void pull(double[] acc,double t) 
	{
		for(int i=0;i<dim;i++)
		{
			v[i]+=t*acc[i];
			loc[i]+=t*v[i];
		}	
	}
	
	//To be used to merge two objects into one(when too close...)
	public boolean merge(Particle object,Cluster cluster, boolean zro) 
	{	//System.out.print("inertia:"+inertia[0]*mass+"+"+object.inertia[0]*object.mass);
		boolean out=false;
		double totalmass=mass+object.mass,absmass=Math.abs(mass)+Math.abs(object.mass);
		charge+=object.charge;
		radius=Math.pow(Math.pow(radius, dim)+Math.pow(object.radius, dim), 1.0/dim);
		if(zro&&totalmass==0)
		{
			for(int i=0;i<dim;i++)
			{
			loc[i]=loc[i]+object.loc[i];
			loc[i]/=2;
			}
			color=Cluster.blur(color,object.color,0.5);
			out=true;
		}
		else if(zro)
		{
			for(int i=0;i<dim;i++)
			{
				loc[i]=object.loc[i];
				v[i]=v[i]+object.v[i]*(object.mass);
				v[i]/=totalmass;
			}
			
			color=Cluster.blur(color,object.color,Math.abs(mass/totalmass));
			
			mass=totalmass;
		}
		else if(totalmass!=0)
		{	
		for(int i=0;i<dim;i++)
		{
			loc[i]=loc[i]*Math.abs(mass)+object.loc[i]*Math.abs(object.mass);
			loc[i]/=absmass;
			if(cluster.relativistic)v[i]+=object.v[i];
			else {
			v[i]=v[i]*(mass)+object.v[i]*(object.mass);
			v[i]/=totalmass;
			}
		}
		
		color=Cluster.blur(color,object.color,Math.abs(mass)/(Math.abs(mass)+Math.abs(object.mass)));
		
		mass=totalmass;
		}
		else if(mass==0)
		{
			for(int i=0;i<dim;i++)
			{
				loc[i]=loc[i]+object.loc[i];
				loc[i]/=2;
				v[i]=v[i]+object.v[i];
				v[i]/=2;
			}
			
			color=Cluster.blur(color,object.color,0.5);
			
		}
		else 
		{
			for(int i=0;i<dim;i++)
			{
				loc[i]=loc[i]*Math.abs(mass)+object.loc[i]*Math.abs(object.mass);
				loc[i]/=(Math.abs(mass)+Math.abs(object.mass));
				v[i]=v[i]*(mass)+object.v[i]*(object.mass);
			    
			}
			
			color=Cluster.blur(color,object.color,0.5);
			
			mass=totalmass;
			out=true;
		}System.out.print("merge "+name+" with "+object.name);
	
		if(name.equals("")&&object.name.equals("")) {if(Math.abs(mass)>Universe.ten) {System.out.print("new name=");name=randomSyllable();System.out.println(name);cluster.screens.add(new Screen( new BufferedImage(Screen.width,Screen.height,BufferedImage.TYPE_3BYTE_BGR),this));}else System.out.println("mass="+mass);}
		else 
		{
			if(name.equals("")&&!object.name.equals(""))cluster.replaceScreen(object,this);
			
			//we don't want the names to get too long://ToDo
			if(name.length()+object.name.length()>rand.nextDouble()*30)
			{
				//while(name.substring(name.length()-1).)
			}
			
			name=name.concat(object.name);
			System.out.println("="+name);
		}//System.out.println("="+inertia[0]*mass);
		
		return out;
	}
	private String randomSyllable() 
	{	
		final int cons=21;
		final String[] letters= {"p","b","m","f","v","t","d","n","th","dh","k","g","ng","ch","h","s","z","sh","zh","l","r","a","e","i","o","u","ae","oe","y","ai","au","ei","ja","je","jo","ju"};
		final int[][]standard= {{0,1,2,3,4},{5,6,7,8,9},{10,11,12,13,14}},
							hiss= {{15,16},{17,18}};
		final int[]special= {19,20},
							vowel= {21,22,23,24,25,26,27,28,29,30,31,32,33,34,35};
		final double[] startl= {0.25,0.4,0.3,0.05},//probability of i consonants at beginning
				twoconsstart={0.1 /*"pf"*/,0.1 /*"pr"*/, 0.2/*"ps"*/,
			 0.15/*"fr"*/,
			0.2 /*"sp"*/, 0.15  /*"sm"*/,0.1 /*"sr"*/},
				threeconsstart= {0.5/*"spr"*/,0.5 /*"spl"*/},
			 endl= {};
		int[]contains;
		final ArrayList<Integer> containsSt=new ArrayList<Integer>();
		containsSt.add(0);/*,con1,2,3,4,5,6,7,8,9,11,12);if (containsSt.*/
		String out="";
		
		//beginning:
		double start=rand.nextDouble();
		int clr = 0, voice,ssh = 0,lr = 0;
		if(start>=startl[0])
		{
			start-=startl[0];
			if(start>=startl[1])
			{
				start-=startl[1];
				if(start>=startl[2])
				{
					start-=startl[2];
					//3consonants
					{
						start/=startl[3];
						start*=standard.length;
						clr=(int)start;
						start=start%1;
						start*=2;
						voice=(int)start;
						start=start%1;
						start*=2;
						ssh=(int)start;
						start=start%1;
						start*=2;
						lr=(int)start;
						start=start%1;
						out=out.concat(letters[hiss[ssh][voice]]+letters[standard[clr][voice]]+letters[special[lr]]);
					}
				}
				else//two consonants
				{
					boolean done=false;
					int i=0;
					start/=startl[2];
					
					while(!done)
					{
						
					if(start<twoconsstart[i])
					{
						done=true;
						start/=twoconsstart[i];
						
						if(i!=6)
						{
							start*=standard.length;
							clr=(int)start;
							start=start%1;
						}
						start*=2;
						voice=(int)start;
						start=start%1;
						if(i==3||(i>4))
						{
							start*=2;
							ssh=(int)start;
							start=start%1;
						}
						if(i==1||i==3||i==6)
						{
							start*=2;
							lr=(int)start;
							start=start%1;
						}
						if(i<3)
							out.concat(letters[standard[clr][voice]]);
						else if(i==3)
							out=out.concat(letters[standard[clr][voice+3]]);
						else
							out=out.concat(letters[hiss[ssh][voice]]);
						//2ndletter
						if(i==4)
							out=out.concat(letters[standard[clr][voice]]);
						else if(i==5)
							out=out.concat(letters[standard[clr][2]]);
						else if(i==0)
							out=out.concat(letters[standard[clr][voice+3]]);
						else if (i==2)
							out=out.concat(letters[hiss[ssh][voice]]);
						else out=out.concat(letters[special[lr]]);
					}
					else
					{
						start-=twoconsstart[i];
						i++;
					}
					}
						
				}
			}
			else//length1 starts
			{
				start*=cons/startl[1];
				out=out.concat(letters[(int)start]);
				start=start%1;
			}
		}
		//vowel
		start*=vowel.length;
		out=out.concat(letters[vowel[(int)start]]);
		start=start%1;
		//end
		if(start>startl[0])
		{
			start-=startl[0];
			if(start>startl[1])
			{
				start-=startl[1];
				if(start>startl[2])
				{
					start-=startl[2];
					//3consonants
					{
						start/=startl[3];
						start*=standard.length;
						clr=(int)start;
						start=start%1;
						start*=2;
						voice=(int)start;
						start=start%1;
						start*=2;
						ssh=(int)start;
						start=start%1;
						start*=2;
						lr=(int)start;
						start=start%1;
						out=out.concat(letters[special[lr]]+letters[standard[clr][voice]]+letters[hiss[ssh][voice]]);
					}
				}
				else//two consonants
				{
					boolean done=false;
					int i=0;
					start/=startl[2];
					while(!done)
					{
						//System.out.println(i+","+start);
					if(start<twoconsstart[i])
					{
						done=true;
						start/=twoconsstart[i];
						
						if(i!=6)
						{
							start*=standard.length;
							clr=(int)start;
							start=start%1;
						}
						start*=2;
						voice=(int)start;
						start=start%1;
						if(i==3||(i>4))
						{
							start*=2;
							ssh=(int)start;
							start=start%1;
						}
						if(i==1||i==3||i==6)
						{
							start*=2;
							lr=(int)start;
							start=start%1;
						}
						if(i<3)
							out=out.concat(letters[standard[clr][voice]]);
						else if(i==3)
							out=out.concat(letters[standard[clr][voice+3]]);
						else
							out=out.concat(letters[hiss[ssh][voice]]);
						//2ndletter
						if(i==4)
							out=out.concat(letters[standard[clr][voice]]);
						else if(i==5)
							out=out.concat(letters[standard[clr][2]]);
						else if(i==0)
							out=out.concat(letters[standard[clr][voice+3]]);
						else if (i==2)
							out=out.concat(letters[hiss[ssh][voice]]);
						else out=out.concat(letters[special[lr]]);
					}
					else
					{
						start-=twoconsstart[i];
						i++;
					}
					}
						
				}
			}
			else//length1 starts
			{
				start*=cons/startl[1];
				out=out.concat(letters[(int)start]);
			}
		}
		return out;
	}
	void saturize() 
	{
		Color c=new Color(color);
		if (!c.equals(new Color(127,127,127)))
		{
		double red=c.getRed()-127.5,green=c.getGreen()-127.5,blue=c.getBlue()-127.5,max=Math.max(Math.abs(red), Math.max(Math.abs(green), Math.abs(blue)));
		red=red*127.5/max+128;green=green*127.5/max+128;blue=blue*127.5/max+128;
		color=new Color((int)red,(int)green,(int)blue).getRGB();
		}
	}
	@Override
	public int compareTo(java.lang.Object o) 
	{
		Particle ob=(Particle)o;
		return (int)Math.signum(mass-ob.mass);
	}

	public void print() 
	{
		System.out.print("loc{");
		for(int i=0;i<dim;i++)
			System.out.print(loc[i]+", ");
		System.out.print("}, v{");
		for(int i=0;i<dim;i++)
			System.out.print(v[i]+", ");
		System.out.println("}, m="+mass+", r="+radius);
		
	}
	//If the force is a first derivative instead of second.
	public void shove(double[] acc, double t) 
	{
		for(int i=0;i<dim;i++)
		{
			loc[i]+=t*acc[i];
			v[i]=acc[i];
		}	
		
	}
	public double torusdistance(double[]dir, double norm,double[]size)
	{
		double out=0;
		for(int i=0;i<dim;i++)
		{
			out+=Math.pow(Math.abs(dir[i]), norm);
		}
		out=Math.pow(out, 1/norm);
		
		return out;
	}
	public void setCharge(double e)
	{
		charge=e;
	}
	public double[] gravAcc(Particle obj)
	{
		double[]acc=new double[dim];
		double d=distance(obj,Cluster.metric);
		for(int k=0;k<dim;k++)
		acc[k]+=-g*obj.mass*Math.pow(d,-(rulescope+1))*(loc[k]-obj.loc[k]);
		
		return acc;
	}
	public double[] torusgravAcc(Particle obj,double d, double[] dir, boolean relativistic)
	{
		double[]acc=new double[dim];
		for(int k=0;k<dim;k++)
		{
			acc[k]+=-g*obj.mass*Math.pow(d,-(rulescope+1))*(dir[k]);
			if(relativistic)
				acc[k]*=mass;
		}
		return acc;
	}
	public double[] elAcc(Particle obj)
	{
		double[]acc=new double[dim];
		double d=distance(obj,Cluster.metric);
		for(int k=0;k<dim;k++)
		acc[k]+=kc*obj.charge*charge/mass*Math.pow(d,-(rulescope+1))*(loc[k]-obj.loc[k]);
		
		return acc;
	}
	
	public double[] relgravForce(Particle obj)
	{
		double[]force=new double[dim];
		double d=distance(obj,Cluster.metric);
		for(int k=0;k<dim;k++)
		force[k]+=-g*mass*obj.mass*Math.pow(d,-(rulescope+1))*(loc[k]-obj.loc[k]);
		//System.out.println("Force="+force[0]+","+force[1]);
		return force;
	}
	public void relPull(double[] force, double t) 
	{
	
		for(int i=0;i<dim;i++)
		{
			v[i]+=t*force[i];
		}
		double pnorm=vdifference(zro(),metric);
		for(int i=0;i<dim;i++)
		{
			double speed=v[i]*c/Math.pow(Math.pow(pnorm, metric)+Math.pow(mass*c, metric), 1.0/metric);
			loc[i]+=t*speed;
			//System.out.println("speed"+speed);
		}
	}
}
