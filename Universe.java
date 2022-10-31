import java.awt.Color;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Random;

//***************************************************
//author: Non-Euclidean Dreamer
// main method displaying the evolution of a universe
//****************************************************
public class Universe 
{ 
	public static boolean relativistic=false,
					toroidal=true,
					manual=false;
	public static String name="classicneg1_";
	public static int background=new Color(0,0,0).getRGB();
	public static int n=500, //number of objects
					dim=2,
					
					counter=0, //starting time
					T=5000, //iterations 
					width=2560,height=2560,
					scale=100,
					ten=100;//what mass is considered considerable 
	public static double t=0.5, //time step
					blur=0.93,//0:only current position visible, 1: whole trace visible
					c=2,//speed of light
					g=1;//gravitational constant
	
	//starting condition for manually set particles
	public static double[][]loc0={{0,0},{10000,0},{10500,0},{1000,0}},
					inertia0={{0,0},{0,1000000000000000.0},{0,3000000000000000.0},{0,-1000},{1.6,1.6}} ;
	public static double[] mass= {2*Math.pow(10,16),Math.pow(10, 12),Math.pow(10, 13),90000,0},
					radius= {1000,200,200,100},
					torussize= {1000,1000};
	public static int[][] color= {{255,0,0},{0,0,255},{0,255,0},{255,255,255}};
	public static BufferedImage screen=new BufferedImage(width,height,BufferedImage.TYPE_4BYTE_ABGR);

	public static void main(String[] args) 
	{ 
		Cluster.relativistic=relativistic;
		Screen.height=height;
		Random rand=new Random();
		Particle.g=g;Particle.c=c;
		Cluster cluster;
		ArrayList<Screen> screens=new ArrayList<Screen>();
		screens.add(new Screen(Particle.zro()));
		
		if(manual)
		{	ArrayList<Particle> obj=new ArrayList<Particle>();
		for(int i=0;i<n;i++)
		{
			int c=rand.nextInt(8); 
			obj.add(new Particle(loc0[i],inertia0[i],mass[i],radius[i],new Color(255*(c%2),255*(c/2%2),255*(c/4%2)).getRGB()));
			
		}
		
		cluster=new Cluster(obj,screens);
		}
		else
		{
			double vmid=01,vsd=0.3;
			double[]zro= {0,0},highloc= {5000,5000},highv= {5,5};
			Particle mid=new Particle(zro,zro,-1,1,background),high=new Particle(highloc,highv,-1,1,background) ,center=new Particle(zro,zro,1,1,background);
			cluster=Cluster.random(n,dim,mid,high);//	Cluster.bigBang(n,vmid,vsd,name);//
			//random: Particle cloud vs bigBang: Start in the middle with outwards velocity
		}
		
		{
			
		cluster.screens.get(0).focus.name=name;
		cluster.centralize();
		
		//cluster.print();
		for(int time=counter;time<T+counter;time++)//for(int k=0;k<7;k++)
		{			
			System.out.println(time+", n="+cluster.objects.size()+", x="+cluster.objects.get(0).loc[0]+", y="+cluster.objects.get(0).loc[1]);
			cluster.torusupdate(t,torussize);
			cluster.torus(torussize); 
			//cluster.centralize();
			//System.out.println("draw"+cluster.screens.size());
			cluster.draw(time);
			
		}
		}
	}

}
