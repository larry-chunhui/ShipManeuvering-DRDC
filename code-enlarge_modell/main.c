#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>


char partName[80];//文件头，

struct Normal
{
	float i;
	float j;
	float k;
};

struct Vertrex
{
	float x;
	float y;
	float z;
};

struct Face
{
	struct Normal normal;
	struct	Vertrex vertex1;
	struct	Vertrex vertex2;
	struct	Vertrex vertex3;
};

//接收三角片面的链表结点
struct Faces
{
	struct	Face data;
	struct	Faces *pointer;
};

int main(int argc,char **argv)
{
	FILE *fp;
	char fileName[64];
	struct Faces *faces = (struct Faces*)malloc(sizeof(struct Faces));//链表的表头
	struct Faces *p = faces;
	struct Faces *q;
	float ratio;
	ratio=0.001;
//	printf("请输入文件名：\n");
//	gets_s(fileName);  //fileName=scanf from the keyboard
//	fileName="tmodel.stl";
	fp = fopen("DTM.stl", "r");


	if (fp)
	{
		char str_a[24];
		char str_b[24];
		char str[24];
		fscanf(fp, "%s", str_a);//该代码的作用仅是移动文件指针，使指针跳过标识符“solid”，将文件指针移到文件名称的位置
		fscanf(fp, "%s", str_a);
		fscanf(fp, "%s", str_a);
		fscanf(fp, "%s", str_a);
		// printf("%s\n", str_a);
		//fscanf(fp, "%s", partName);//读入文件名称,first line in the file
		// printf("%s\n", partName);
		//循环读取三角片面的法线和顶点数据
		while (!feof(fp))
		{
			p->pointer= (struct Faces*)malloc(sizeof(struct Faces));//创建新结点，用于接收数据
			q = p;
			p = p->pointer;
			fscanf(fp, "%s%s\n", str_a,str_b);//使文件指针跳过标识符“facetnormal”，指到法线数据部分
		        // printf("%s %s\n", str_a,str_b);
			fscanf(fp, "%f%f%f", &p->data.normal.i, &p->data.normal.j, &p->data.normal.k);//读取法线数据
			// printf("%f,%f,%f\n", p->data.normal.i, p->data.normal.j, p->data.normal.k);
			//return(0);

			fscanf(fp, "%s%s\n", str_a,str_b);//使文件指针跳过标识符“outer loop”
			fscanf(fp, "%s", str);//使文件指针跳过标识符“vertex”，指到顶点1的数据部分
			fscanf(fp, "%f%f%f", &p->data.vertex1.x, &p->data.vertex1.y, &p->data.vertex1.z);

			fscanf(fp, "%s", str);//使文件指针跳过标识符“vertex”，指到顶点2的数据部分
			fscanf(fp, "%f%f%f", &p->data.vertex2.x, &p->data.vertex2.y, &p->data.vertex2.z);
			fscanf(fp, "%s", str);//使文件指针跳过标识符“vertex”，指到顶点3的数据部分
			fscanf(fp, "%f%f%f", &p->data.vertex3.x, &p->data.vertex3.y, &p->data.vertex3.z);
			fscanf(fp, "%s", str);//使文件指针跳过标识符“endloop”
			fscanf(fp, "%s", str);//使文件指针跳过标识符“endfacet“
		}
		free(q->pointer);//由于文件末尾有一行字符”endsolid .....“，造成多创建了一个结点，该结点的数据是无用的，因此要释放
		q->pointer = NULL;//最后一个结点的数据释放完之后，将指向最后结点的指针置为NULL，方便遍历链表

	}
	else
	{
		printf("打开文件失败！");
	}
	fclose(fp);
// output file
	p = faces->pointer;//获取链表第二个结点的指针
	fp=fopen("DTM-scaled.stl","w");
	fprintf(fp,"%s\n","solid Created by Gmsh");
	while (p!=NULL)
	{
		//循环输出链表
//		printf("法线：\n");
		fprintf(fp,"%s %f  %f  %f\n","facet normal", p->data.normal.i, p->data.normal.j, p->data.normal.k);
 		fprintf(fp,"%s\n","outer loop");
		fprintf(fp,"%s %f  %f  %f\n","	vertex", p->data.vertex1.x*ratio, p->data.vertex1.y*ratio, p->data.vertex1.z*ratio);
		fprintf(fp,"%s %f  %f  %f\n","	vertex", p->data.vertex2.x*ratio, p->data.vertex2.y*ratio, p->data.vertex2.z*ratio);
		fprintf(fp,"%s %f  %f  %f\n","	vertex", p->data.vertex3.x*ratio, p->data.vertex3.y*ratio, p->data.vertex3.z*ratio);
 		fprintf(fp,"%s\n","endloop");
		fprintf(fp,"%s\n","endfacet");

		p = p->pointer;
	}
		fprintf(fp,"%s\n","endsolid Created by Gmsh");
	fclose(fp);
	return 0;
//--------------------------------------------------------
	p = faces->pointer;//获取链表第二个结点的指针
	while (p!=NULL)
	{
		//循环输出链表
		printf("法线：\n");
		printf("%f  %f  %f\n", p->data.normal.i, p->data.normal.j, p->data.normal.k);
		printf("顶点：\n");
		printf("%f  %f  %f\n", p->data.vertex1.x*ratio, p->data.vertex1.y*ratio, p->data.vertex1.z*ratio);
		printf("%f  %f  %f\n", p->data.vertex2.x*ratio, p->data.vertex2.y*ratio, p->data.vertex2.z*ratio);
		printf("%f  %f  %f\n", p->data.vertex3.x*ratio, p->data.vertex3.y*ratio, p->data.vertex3.z*ratio);
		p = p->pointer;
	}

	return 0;

}












