import numpy as np
import matplotlib.pyplot as plt

strain_HorCyl_0, stress_HorCyl_0, strain_HorCyl_1, stress_HorCyl_1, strain_HorCyl_2, stress_HorCyl_2 = \
    np.loadtxt('HorCylinder_data.csv', delimiter=',', unpack=True)

strain_VertCyl_0, stress_VertCyl_0, strain_VertCyl_1, stress_VertCyl_1, strain_VertCyl_2, stress_VertCyl_2 = \
    np.loadtxt('VertCylinder_data.csv', delimiter=',', unpack=True)

strain_Sphere_0, stress_Sphere_0, strain_Sphere_1, stress_Sphere_1, strain_Sphere_2, stress_Sphere_2, strain_Sphere_3,\
            stress_Sphere_3 = \
    np.loadtxt('Sphere_data.csv', delimiter=',', unpack=True)

strain_Plates_0, stress_Plates_0, strain_Plates_1, stress_Plates_1, strain_Plates_2, stress_Plates_2, strain_Plates_3,\
            stress_Plates_3 = \
    np.loadtxt('Plates_data.csv', delimiter=',', unpack=True)


print('Youngs Modulus for pure matrix = ', stress_HorCyl_0[-1]/strain_HorCyl_0[-1]/1e09, 'GPa')
print('Youngs Modulus for horizontal cylinders occupying 1% volume fraction = ', stress_HorCyl_1[-1]/strain_HorCyl_1[-1]/1e09, 'GPa')
print('Youngs Modulus for horizontal cylinders occupying 2% volume fraction = ', stress_HorCyl_2[-1]/strain_HorCyl_2[-1]/1e09, 'GPa')
print('Youngs Modulus for vertical cylinders occupying 1% volume fraction = ', stress_VertCyl_1[-1]/strain_VertCyl_1[-1]/1e09, 'GPa')
print('Youngs Modulus for vertical cylinders occupying 2% volume fraction = ', stress_VertCyl_2[-1]/strain_VertCyl_2[-1]/1e09, 'GPa')
print('Youngs Modulus for spheres occupying 1% volume fraction = ', stress_Sphere_1[-1]/strain_Sphere_1[-1]/1e09, 'GPa')
print('Youngs Modulus for spheres occupying 2% volume fraction = ', stress_Sphere_2[-1]/strain_Sphere_2[-1]/1e09, 'GPa')
print('Youngs Modulus for spheres occupying 3% volume fraction = ', stress_Sphere_3[-1]/strain_Sphere_3[-1]/1e09, 'GPa')
print('Youngs Modulus for plates occupying 1% volume fraction = ', stress_Plates_1[-1]/strain_Plates_1[-1]/1e09, 'GPa')
print('Youngs Modulus for plates occupying 2% volume fraction = ', stress_Plates_2[-1]/strain_Plates_2[-1]/1e09, 'GPa')
print('Youngs Modulus for plates occupying 3% volume fraction = ', stress_Plates_3[-1]/strain_Plates_3[-1]/1e09, 'GPa')
plt.plot(strain_HorCyl_0[1:], stress_HorCyl_0[1:]/1e09, 'r')
plt.plot(strain_HorCyl_1[1:], stress_HorCyl_1[1:]/1e09, 'g')
plt.plot(strain_HorCyl_2[1:], stress_HorCyl_2[1:]/1e09, 'b')
plt.xlabel('True Strain')
plt.ylabel('True Stress (in GPa)')
plt.legend(('Pure solid', '1% composite', '2% composite'))
plt.title('Cylinders arranged perpendicular to loading direction')
plt.show()

plt.plot(strain_VertCyl_0[1:], stress_VertCyl_0[1:]/1e09, 'r')
plt.plot(strain_VertCyl_1[1:], stress_VertCyl_1[1:]/1e09, 'g')
plt.plot(strain_VertCyl_2[1:], stress_VertCyl_2[1:]/1e09, 'b')
plt.xlabel('True Strain')
plt.ylabel('True Stress (in GPa)')
plt.legend(('Pure solid', '1% composite', '2% composite'))
plt.title('Cylinders arranged along the loading direction')
plt.show()

plt.plot(strain_Sphere_0[1:], stress_Sphere_0[1:]/1e09, 'r')
plt.plot(strain_Sphere_1[1:], stress_Sphere_1[1:]/1e09, 'g')
plt.plot(strain_Sphere_2[1:], stress_Sphere_2[1:]/1e09, 'b')
plt.plot(strain_Sphere_3[1:], stress_Sphere_3[1:]/1e09, 'brown')
plt.xlabel('True Strain')
plt.ylabel('True Stress (in GPa)')
plt.legend(('Pure solid', '1% composite', '2% composite', '3% composite'))
plt.title('Spherical Inclusions')
plt.show()

plt.plot(strain_Plates_0[1:], stress_Plates_0[1:]/1e09, 'r')
plt.plot(strain_Plates_1[1:], stress_Plates_1[1:]/1e09, 'g')
plt.plot(strain_Plates_2[1:], stress_Plates_2[1:]/1e09, 'b')
plt.plot(strain_Plates_3[1:], stress_Plates_3[1:]/1e09, 'brown')
plt.xlabel('True Strain')
plt.ylabel('True Stress (in GPa)')
plt.legend(('Pure solid', '1% composite', '2% composite', '3% composite'))
plt.title('Plate-like Inclusions')
plt.show()

plt.plot(strain_HorCyl_1[1:], stress_HorCyl_1[1:]/1e09)
plt.plot(strain_VertCyl_1[1:], stress_VertCyl_1[1:]/1e09)
plt.plot(strain_Sphere_1[1:], stress_Sphere_1[1:]/1e09)
plt.plot(strain_Plates_1[1:], stress_Plates_1[1:]/1e09)
plt.xlabel('True Strain')
plt.ylabel('True Stress (in GPa)')
plt.legend(('Horizontal cylinders', 'Vertical cylinders', 'Spheres', 'Plates'))
plt.title('Different types of inclusions at 1% volume fraction')
plt.show()

plt.plot(strain_HorCyl_2[1:], stress_HorCyl_2[1:]/1e09)
plt.plot(strain_VertCyl_2[1:], stress_VertCyl_2[1:]/1e09)
plt.plot(strain_Sphere_2[1:], stress_Sphere_2[1:]/1e09)
plt.plot(strain_Plates_2[1:], stress_Plates_2[1:]/1e09)
plt.xlabel('True Strain')
plt.ylabel('True Stress (in GPa)')
plt.legend(('Horizontal cylinders', 'Vertical cylinders', 'Spheres', 'Plates'))
plt.title('Different types of inclusions at 2% volume fraction')
plt.show()

plt.plot(strain_Sphere_3[1:], stress_Sphere_3[1:]/1e09)
plt.plot(strain_Plates_3[1:], stress_Plates_3[1:]/1e09)
plt.xlabel('True Strain')
plt.ylabel('True Stress (in GPa)')
plt.legend(('Spheres', 'Plates'))
plt.title('Different types of inclusions at 3% volume fraction')
plt.show()