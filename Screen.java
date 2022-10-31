import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

//**********************************************************************************************
//Author: Non-Euclidean Dreamer
//A Screen displaying a Cluster
//**********************************************************************************************


public class Screen 
{
	static int width=2560, height=2560, background=Color.black.getRGB();
	static double blur=0.99 ,metric=Cluster.metric;
	 double scale=100;
	BufferedImage scr;
	Particle focus;
	
	public Screen(Particle f)
	{
		scr=new BufferedImage(width,height, BufferedImage.TYPE_3BYTE_BGR);
		for(int i=0;i<width;i++)
			for(int j=0;j<height;j++)
				scr.setRGB(i,j, background);
		focus=f;
		scale=Math.max(f.radius*100,500);
	}
	
	public Screen(BufferedImage s,Particle f)
	{
		scr=s;
		focus=f;
		scale=Math.max(f.radius*100,500);
	}

	public void draw(Cluster cluster,int t)
	{
		for(int i=0;i<width;i++)
		{
			for(int j=0;j<height;j++)
			{
				scr.setRGB(i, j, Cluster.blur(scr.getRGB(i, j),background,blur));
			}
		}
		//System.out.println("blablub");
		for(int k=0;k<cluster.objects.size();k++)
		{
			Particle obj= cluster.objects.get(k);//System.out.println("x="+obj.loc[0]);
			double[]center= {width/2+(obj.loc[0]-focus.loc[0])/scale*width,height/2+(obj.loc[1]-focus.loc[1])/scale*width};//System.out.println("center"+center[0]+","+center[1]);
			double r=obj.radius*width/scale;//System.out.println("radius="+r);
			for(int i=Math.max((int)(center[0]-r),0);i<center[0]+r+1&&i<width;i++)
			{
				//System.out.println("i="+i);
				double s=Math.pow(Math.pow(r, metric)-Math.pow(Math.abs(center[0]-i), metric), 1.0/metric);
				for(int j=Math.max(0,(int)(center[1]-s));j<center[1]+s+1&&j<height;j++)
				{
					//System.out.println(i+","+j+","+obj.color);//here+
					try{scr.setRGB(i,j,obj.color);}catch(ArrayIndexOutOfBoundsException e) {}
				}
			}
				
		} 
		File outputfile = new File(focus.name+t+".png");
		try 
		{
			ImageIO.write(scr, "png", outputfile);
		} catch (IOException e) {
			System.out.println("IOException");
			e.printStackTrace(); 
		}	
	}
}
