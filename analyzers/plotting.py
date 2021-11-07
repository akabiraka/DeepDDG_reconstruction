import matplotlib.pyplot as plt
import numpy as np
from numpy.core.fromnumeric import size

def plot_train_losses(train_losses):    
    plt.plot(train_losses)
    plt.savefig("outputs/images/model_analysis/{}.pdf".format("train_loss"), dpi=300, format="pdf", bbox_inches='tight', pad_inches=0.0)
    plt.cla()
    # plt.show()

def plot_pred_vs_exp_ddg(exp_ddgs, pred_ddgs):
    x = np.array(exp_ddgs)
    y = np.array(pred_ddgs)
    plt.scatter(x, y, c="salmon", marker=".")

    m, b = np.polyfit(x, y, 1) #m = slope, b=intercept
    plt.plot(x, m*x + b, c="lightgreen")

    plt.xlabel("Expected ddG (kcal/mol)")
    plt.ylabel("Predicted ddG (kcal/mol)")
    plt.savefig("outputs/images/model_analysis/{}.pdf".format("pred_vs_exp_ddg"), dpi=300, format="pdf", bbox_inches='tight', pad_inches=0.0)
    plt.cla()
    # plt.show()
    
def plot_deviation_from_expected(exp_ddgs, pred_ddgs):
    x = np.array(exp_ddgs)
    y = np.array(pred_ddgs)
    error = x-y
    
    zeros = np.zeros_like(x)
    plt.errorbar(x, zeros, yerr=error, fmt='o', color='green',
             ecolor='salmon', elinewidth=3, capsize=0, alpha=.5)

    plt.xlabel("Expected ddG (kcal/mol)")
    plt.ylabel("Error deviation (kcal/mol)")
    plt.savefig("outputs/images/model_analysis/{}.pdf".format("deviation_from_expected"), dpi=300, format="pdf", bbox_inches='tight', pad_inches=0.0)
    plt.cla()
    # plt.show()
        
train_losses = [0.0179902296513319, 0.016692524775862694, 0.015788046643137932, 0.014784918166697025, 0.013171518221497536, 0.012993525713682175, 0.011737145483493805, 0.00948218535631895, 0.008877377025783062, 0.007874509319663048]
exp_ddgs= [-2.549999952316284, -3.5, -1.899999976158142, -3.049999952316284, -0.699999988079071, -2.0, -0.10000000149011612, 0.6600000262260437, 0.7099999785423279, 0.0, 0.17000000178813934, 0.30000001192092896, -0.30000001192092896, -1.2000000476837158, -2.690000057220459, -1.9800000190734863, -1.7000000476837158, -1.3700000047683716, -2.049999952316284, -3.009999990463257, -0.6000000238418579, -2.880000114440918, -1.3600000143051147, -1.8600000143051147, -2.0199999809265137, -2.880000114440918, -1.9700000286102295, 0.3449999988079071, 0.7655502557754517, 0.5502392053604126, -1.9138755798339844, -1.4832535982131958, -1.6746411323547363, 0.5263158082962036, 0.2631579041481018, 1.6267942190170288, -0.3349282443523407, -0.16099999845027924, -0.2199999988079071, -0.6700000166893005, -1.628342866897583, -1.2065000534057617, -0.15000000596046448, -0.3100000023841858, -2.03507924079895, -1.0140000581741333, -0.29477372765541077, -0.14000000059604645, -0.07999999821186066, -1.5057361125946045, -2.4600000381469727, -0.4461440443992615, 0.8126195073127747, -0.33000001311302185, -2.130000114440918, -1.440000057220459, -4.619999885559082, -0.17000000178813934, -0.5, -0.25999999046325684, 0.23900572955608368, -2.5999999046325684, -1.7999999523162842, -1.100000023841858, -1.100000023841858, -1.100000023841858, -0.20000000298023224, -1.2000000476837158, -0.10000000149011612, -1.2999999523162842, 0.699999988079071, -0.30000001192092896, 0.10000000149011612, -1.100000023841858, -3.700000047683716, -1.2999999523162842, -2.799999952316284, 0.5, 0.5, -0.5, -1.399999976158142, -2.299999952316284, -3.6500000953674316, -1.600000023841858, -2.799999952316284, -0.6000000238418579, -2.0, -0.6000000238418579, 0.0, 0.4000000059604645, -0.5, -0.10000000149011612, -1.399999976158142, -0.4000000059604645, -1.2999999523162842, -0.20000000298023224, -2.799999952316284, -3.799999952316284, -2.4000000953674316, -2.799999952316284, -3.799999952316284, -2.200000047683716, 0.20000000298023224, -0.4000000059604645, 0.4000000059604645, 0.4000000059604645, -0.699999988079071, 0.5, -1.2999999523162842, -0.699999988079071, -1.7999999523162842, 0.0, -1.2999999523162842, 0.0, -0.30000001192092896, -1.899999976158142, 0.30000001192092896, -0.6000000238418579, 0.5, -0.5, 0.30000001192092896, -0.10000000149011612, -0.20000000298023224, -0.4000000059604645, -2.299999952316284, -1.1100000143051147, -2.0999999046325684, 1.0, -3.5, 0.20000000298023224, -1.0, -1.600000023841858, -5.0, -2.200000047683716, -1.600000023841858, -1.7999999523162842, -3.9000000953674316, -2.299999952316284, -1.600000023841858, -2.5, 3.200000047683716, -1.2999999523162842, -1.399999976158142, -0.800000011920929, -2.4000000953674316, -0.20000000298023224, -1.6507177352905273, -3.4688994884490967, -1.5071769952774048, -1.2918660640716553, -4.330143451690674, -1.7224880456924438, -1.4114832878112793, -1.7703349590301514, -5.358851909637451, -1.6746411323547363, -3.9712917804718018, -0.4545454680919647, -1.6028708219528198, -1.8421052694320679, -1.4114832878112793, -2.6315789222717285, -1.8421052694320679, -0.2870813310146332, -2.5358850955963135, -1.6267942190170288, -1.8660286664962769, -1.8899521827697754, -0.023923445492982864, -1.5071769952774048, 0.43062201142311096, -2.2248804569244385, -0.23923444747924805, -0.800000011920929, -2.7899999618530273, -0.5, -1.190000057220459, 0.4000000059604645, -0.75, -0.5, -1.690000057220459, -0.30000001192092896, -2.4000000953674316, -2.5899999141693115, -2.700000047683716, -2.2249999046325684, 0.30000001192092896, -2.0, 0.20000000298023224, -1.899999976158142, -0.6000000238418579, -2.450000047683716, -1.2000000476837158, 0.10000000149011612, 0.0, -1.600000023841858, -1.7999999523162842, -1.159999966621399, 0.4000000059604645, -2.700000047683716, -1.725000023841858, -1.7999999523162842, -3.4000000953674316]
pred_ddgs= [-1.9486194849014282, -2.4017463624477386, -1.4583303034305573, -1.1819985508918762, -1.299278885126114, -1.340779960155487, -0.6423831731081009, 0.7766865193843842, 0.3544605150818825, -0.21026849746704102, 0.18652979284524918, 0.3328954428434372, 0.33526964485645294, -0.967174619436264, -1.3198046386241913, -1.6654270887374878, -1.541711688041687, -1.1885006725788116, -1.3435500860214233, -1.550724357366562, -0.8000461012125015, -1.6917979717254639, -0.5867435783147812, -0.9396127611398697, -1.1140385270118713, -1.6898338496685028, -0.9858979284763336, -0.9321232885122299, 0.44287804514169693, 0.4328112304210663, -1.3939681649208069, -1.5383581817150116, -1.5353992581367493, 0.6122950464487076, 0.6081249192357063, 1.1340329051017761, -0.2683468349277973, 0.03062037518247962, 0.027354282792657614, -0.8679684996604919, -1.152164712548256, -0.8774654567241669, -0.8758073300123215, 0.2015729434788227, -1.4567261934280396, -0.3023291565477848, 0.18751714378595352, -0.19159551709890366, -0.13625339604914188, -1.7109471559524536, -1.6970933973789215, -0.31193677335977554, -0.29427655041217804, -0.5055642127990723, -1.6378980875015259, -0.7481106370687485, -2.4442386627197266, -0.4186943918466568, -0.7420635968446732, -0.16873011365532875, 0.2349518984556198, -1.560661494731903, -1.3000604510307312, -0.8620009571313858, -0.5928007885813713, -0.9237552434206009, -0.7029232382774353, -0.7037480920553207, -0.28200818225741386, -0.2837328240275383, -0.08451823145151138, -0.0864078477025032, -0.23065727204084396, -0.2230655774474144, -1.9882604479789734, -1.9758939743041992, -0.9451572597026825, -0.936841294169426, 0.25961291044950485, 0.25691457092761993, -0.7685965299606323, -1.6662967205047607, -1.6152320802211761, -1.6154088079929352, -1.5647777915000916, -1.5634678304195404, -1.2045644223690033, -1.1989396065473557, -0.9532289206981659, -0.14651041477918625, -0.14475835487246513, -0.822489783167839, -0.8249910920858383, -0.6331422924995422, -0.6311920285224915, 0.16656268388032913, -1.550731509923935, -1.6925734281539917, -1.6920189559459686, -1.9160951673984528, -1.6429969668388367, -2.0196932554244995, -0.183821190148592, -0.18279772251844406, -0.7902836054563522, 0.2360539697110653, 0.2370777539908886, -0.037630125880241394, -0.0366451614536345, -0.680161789059639, -0.6834220886230469, -0.9897654503583908, -0.9880901873111725, 0.07062749937176704, -1.1902394145727158, -1.1897853761911392, 0.18063675612211227, 0.18220704048871994, 0.24401290342211723, 0.2429644577205181, 0.7239982485771179, 0.7253150641918182, -1.2935999035835266, -0.22065360099077225, -1.7933200299739838, -0.6806939095258713, -1.58599391579628, 0.16183599829673767, -1.1165518313646317, -1.1028394848108292, -1.1205150187015533, -1.5561926364898682, -2.041127383708954, -1.9791589677333832, -1.8819805979728699, -1.730717122554779, -2.0545952022075653, -1.9941845536231995, -1.5803661942481995, -2.2519367933273315, 2.31476753950119, -1.1358537524938583, -1.1854398250579834, -1.3023664057254791, -1.4074039459228516, 0.01748063717968762, -0.8639592677354813, -2.4059809744358063, -1.0415199398994446, -0.7323930412530899, -3.130461871623993, -1.1309249699115753, -1.0782834887504578, -1.4252348244190216, -3.6378678679466248, -1.0712291300296783, -1.8315055966377258, -0.20206114277243614, -1.1644978821277618, -1.9459405541419983, -1.2444230169057846, -2.1725504100322723, -0.7611279934644699, -0.5352237820625305, -1.9020962715148926, -1.4049673080444336, -1.4058278501033783, -1.1645717918872833, 0.7843435555696487, -0.8233178406953812, -0.8308973908424377, -0.8468014001846313, -0.9170816838741302, -1.4343653619289398, -1.4339147508144379, -0.9771207720041275, -0.9775014221668243, -0.16392672434449196, -0.1653524488210678, -1.5304136276245117, -1.5319626033306122, -0.5595194548368454, -1.8611979484558105, -1.8569053709506989, -1.315131038427353, -1.50190070271492, -0.8152244240045547, -1.449878215789795, -0.6544753164052963, -1.4031164348125458, -1.2379127740859985, -1.239597499370575, -1.2375861406326294, 0.18541550263762474, -0.8833858370780945, -0.8825539797544479, -0.5116759613156319, -0.5134740471839905, 0.37921641021966934, -1.3684676587581635, -1.0493389517068863, -1.0521796345710754, -2.29976624250412]
    
# plot_train_losses(train_losses)
# plot_pred_vs_exp_ddg(exp_ddgs, pred_ddgs)
plot_deviation_from_expected(exp_ddgs, pred_ddgs)