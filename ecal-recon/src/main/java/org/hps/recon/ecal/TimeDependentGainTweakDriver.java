package org.hps.recon.ecal;

import java.util.ArrayList;
//import java.lang.reflect.Modifier;;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import org.hps.rundb.RunManager;
import org.lcsim.event.CalorimeterHit;
import org.lcsim.event.Cluster;
import org.lcsim.event.EventHeader;
import org.lcsim.event.ReconstructedParticle;
import org.lcsim.event.Vertex;
import org.lcsim.event.base.BaseCluster;
import org.lcsim.event.base.BaseReconstructedParticle;
import org.lcsim.event.base.BaseVertex;
import org.lcsim.util.Driver;

import hep.physics.vec.BasicHepLorentzVector;
import hep.physics.vec.Hep3Vector;
/**
 * Driver for tweaking precooked data with time dependent gains on the ecal.  
 * @author spaul
 *
 */
public abstract class TimeDependentGainTweakDriver extends Driver {



    
    private String ecalClusterCollection = "EcalClustersCorr";
    //temporary patch for a bug that affects the timestamps in the files.  
    long ts_start[] = new long[8101];
    TimeDependentGainTweakDriver(){
        ts_start[5254] = 1430531592;
        ts_start[5255] = 1430532260;
        ts_start[5256] = 1430533541;
        ts_start[5257] = 1430535151;
        ts_start[5258] = 1430542206;
        ts_start[5259] = 1430542675;
        ts_start[5263] = 1430546354;
        ts_start[5264] = 1430547410;
        ts_start[5265] = 1430548544;
        ts_start[5266] = 1430552110;
        ts_start[5267] = 1430552721;
        ts_start[5268] = 1430554181;
        ts_start[5269] = 1430554640;
        ts_start[5270] = 1430555077;
        ts_start[5271] = 1430557345;
        ts_start[5272] = 1430558496;
        ts_start[5273] = 1430559758;
        ts_start[5274] = 1430561958;
        ts_start[5275] = 1430563735;
        ts_start[5278] = 1430569003;
        ts_start[5279] = 1430573659;
        ts_start[5280] = 1430575235;
        ts_start[5281] = 1430576624;
        ts_start[5282] = 1430577077;
        ts_start[5283] = 1430578077;
        ts_start[5284] = 1430579155;
        ts_start[5285] = 1430580057;
        ts_start[5286] = 1430580574;
        ts_start[5287] = 1430584864;
        ts_start[5288] = 1430585487;
        ts_start[5292] = 1430589041;
        ts_start[5293] = 1430593882;
        ts_start[5294] = 1430595169;
        ts_start[5295] = 1430597419;
        ts_start[5296] = 1430598236;
        ts_start[5297] = 1430599368;
        ts_start[5298] = 1430600685;
        ts_start[5299] = 1430603416;
        ts_start[5302] = 1430606506;
        ts_start[5303] = 1430610042;
        ts_start[5304] = 1430612340;
        ts_start[5305] = 1430613347;
        ts_start[5306] = 1430617294;
        ts_start[5307] = 1430622752;
        ts_start[5309] = 1430625710;
        ts_start[5310] = 1430626007;
        ts_start[5311] = 1430631240;
        ts_start[5312] = 1430634332;
        ts_start[5313] = 1430637642;
        ts_start[5314] = 1430643044;
        ts_start[5315] = 1430648262;
        ts_start[5316] = 1430653144;
        ts_start[5318] = 1430657483;
        ts_start[5319] = 1430657999;
        ts_start[5322] = 1430660898;
        ts_start[5329] = 1430667372;
        ts_start[5330] = 1430669669;
        ts_start[5331] = 1430674712;
        ts_start[5332] = 1430678518;
        ts_start[5333] = 1430679299;
        ts_start[5334] = 1430681900;
        ts_start[5339] = 1430686134;
        ts_start[5340] = 1430690859;
        ts_start[5341] = 1430692293;
        ts_start[5342] = 1430693251;
        ts_start[5344] = 1430695366;
        ts_start[5345] = 1430702311;
        ts_start[5346] = 1430706553;
        ts_start[5347] = 1430707919;
        ts_start[5348] = 1430717913;
        ts_start[5349] = 1430721379;
        ts_start[5350] = 1430723626;
        ts_start[5351] = 1430729008;
        ts_start[5375] = 1430792995;
        ts_start[5378] = 1430795511;
        ts_start[5379] = 1430796804;
        ts_start[5380] = 1430797414;
        ts_start[5381] = 1430802383;
        ts_start[5382] = 1430807414;
        ts_start[5383] = 1430808710;
        ts_start[5384] = 1430810393;
        ts_start[5385] = 1430815187;
        ts_start[5386] = 1430816114;
        ts_start[5387] = 1430819470;
        ts_start[5402] = 1430869980;
        ts_start[5403] = 1430874952;
        ts_start[5404] = 1430875594;
        ts_start[5405] = 1430882235;
        ts_start[5406] = 1430885395;
        ts_start[5407] = 1430887782;
        ts_start[5409] = 1430890541;
        ts_start[5410] = 1430890843;
        ts_start[5411] = 1430896120;
        ts_start[5412] = 1430903027;
        ts_start[5538] = 1431240789;
        ts_start[5541] = 1431242995;
        ts_start[5542] = 1431243817;
        ts_start[5546] = 1431246831;
        ts_start[5547] = 1431247555;
        ts_start[5548] = 1431247963;
        ts_start[5549] = 1431255679;
        ts_start[5550] = 1431263016;
        ts_start[5554] = 1431266372;
        ts_start[5558] = 1431270429;
        ts_start[5559] = 1431277631;
        ts_start[5560] = 1431284830;
        ts_start[5562] = 1431291930;
        ts_start[5563] = 1431293690;
        ts_start[5564] = 1431294788;
        ts_start[5565] = 1431295273;
        ts_start[5566] = 1431295865;
        ts_start[5567] = 1431296903;
        ts_start[5568] = 1431298077;
        ts_start[5569] = 1431305422;
        ts_start[5575] = 1431312286;
        ts_start[5576] = 1431318777;
        ts_start[5577] = 1431321772;
        ts_start[5578] = 1431328812;
        ts_start[5579] = 1431335968;
        ts_start[5597] = 1431409254;
        ts_start[5598] = 1431409599;
        ts_start[5601] = 1431410986;
        ts_start[5602] = 1431411496;
        ts_start[5603] = 1431411845;
        ts_start[5604] = 1431412208;
        ts_start[5605] = 1431412477;
        ts_start[5610] = 1431417220;
        ts_start[5611] = 1431420077;
        ts_start[5615] = 1431475209;
        ts_start[5616] = 1431476398;
        ts_start[5617] = 1431476909;
        ts_start[5618] = 1431477882;
        ts_start[5619] = 1431479565;
        ts_start[5620] = 1431480126;
        ts_start[5621] = 1431480619;
        ts_start[5622] = 1431482260;
        ts_start[5623] = 1431482527;
        ts_start[5624] = 1431487449;
        ts_start[5626] = 1431491664;
        ts_start[5631] = 1431554773;
        ts_start[5632] = 1431555528;
        ts_start[5633] = 1431557310;
        ts_start[5634] = 1431557707;
        ts_start[5635] = 1431563447;
        ts_start[5636] = 1431564402;
        ts_start[5637] = 1431566280;
        ts_start[5638] = 0;
        ts_start[5639] = 1431568924;
        ts_start[5640] = 1431569393;
        ts_start[5641] = 1431570282;
        ts_start[5642] = 1431570500;
        ts_start[5643] = 1431574631;
        ts_start[5644] = 1431578548;
        ts_start[5645] = 1431583298;
        ts_start[5646] = 1431583569;
        ts_start[5648] = 1431584989;
        ts_start[5649] = 1431586081;
        ts_start[5650] = 1431587742;
        ts_start[5651] = 1431589301;
        ts_start[5652] = 1431590647;
        ts_start[5653] = 1431592142;
        ts_start[5654] = 1431593135;
        ts_start[5655] = 1431594701;
        ts_start[5656] = 1431595392;
        ts_start[5657] = 1431597375;
        ts_start[5683] = 0;
        ts_start[5685] = 1431652250;
        ts_start[5686] = 1431652654;
        ts_start[5687] = 0;
        ts_start[5688] = 1431653987;
        ts_start[5689] = 1431654232;
        ts_start[5690] = 0;
        ts_start[5691] = 1431662772;
        ts_start[5692] = 1431664090;
        ts_start[5693] = 1431671748;
        ts_start[5694] = 1431679041;
        ts_start[5695] = 1431684740;
        ts_start[5696] = 1431694000;
        ts_start[5697] = 1431699980;
        ts_start[5698] = 1431701983;
        ts_start[5702] = 0;
        ts_start[5703] = 1431712837;
        ts_start[5704] = 1431713140;
        ts_start[5705] = 1431713956;
        ts_start[5706] = 1431714831;
        ts_start[5707] = 1431719299;
        ts_start[5708] = 1431719497;
        ts_start[5709] = 1431719601;
        ts_start[5710] = 1431719737;
        ts_start[5711] = 1431725802;
        ts_start[5712] = 0;
        ts_start[5713] = 1431729714;
        ts_start[5714] = 1431738481;
        ts_start[5715] = 1431738722;
        ts_start[5717] = 1431746119;
        ts_start[5718] = 0;
        ts_start[5721] = 1431747612;
        ts_start[5722] = 1431747876;
        ts_start[5723] = 1431749092;
        ts_start[5724] = 1431756289;
        ts_start[5725] = 1431763465;
        ts_start[5726] = 1431769094;
        ts_start[5727] = 0;
        ts_start[5728] = 1431773080;
        ts_start[5729] = 1431775499;
        ts_start[5730] = 1431775833;
        ts_start[5731] = 0;
        ts_start[5733] = 1431776178;
        ts_start[5734] = 0;
        ts_start[5735] = 0;
        ts_start[5736] = 0;
        ts_start[5737] = 1431778802;
        ts_start[5738] = 1431780334;
        ts_start[5739] = 1431780468;
        ts_start[5740] = 1431786552;
        ts_start[5741] = 1431786760;
        ts_start[5742] = 1431793915;
        ts_start[5743] = 1431801798;
        ts_start[5744] = 0;
        ts_start[5745] = 1431806648;
        ts_start[5747] = 1431807613;
        ts_start[5748] = 1431808568;
        ts_start[5749] = 1431809108;
        ts_start[5750] = 0;
        ts_start[5751] = 0;
        ts_start[5752] = 1431810286;
        ts_start[5753] = 1431810807;
        ts_start[5754] = 1431811200;
        ts_start[5755] = 1431811947;
        ts_start[5756] = 1431812449;
        ts_start[5757] = 1431812701;
        ts_start[5758] = 0;
        ts_start[5759] = 0;
        ts_start[5760] = 1431814426;
        ts_start[5761] = 1431814942;
        ts_start[5762] = 1431815403;
        ts_start[5763] = 1431815785;
        ts_start[5764] = 1431816062;
        ts_start[5765] = 1431817384;
        ts_start[5766] = 1431817736;
        ts_start[5767] = 1431825302;
        ts_start[5768] = 1431825654;
        ts_start[5769] = 1431832455;
        ts_start[5770] = 1431839492;
        ts_start[5771] = 1431846393;
        ts_start[5772] = 1431853732;
        ts_start[5773] = 1431861164;
        ts_start[5774] = 1431870016;
        ts_start[5775] = 1431871806;
        ts_start[5776] = 1431876070;
        ts_start[5778] = 1431882602;
        ts_start[5779] = 1431883123;
        ts_start[5780] = 0;
        ts_start[5781] = 1431885610;
        ts_start[5782] = 1431888000;
        ts_start[5783] = 1431896031;
        ts_start[5784] = 1431905102;
        ts_start[5785] = 1431905926;
        ts_start[5786] = 1431906763;
        ts_start[5788] = 1431909065;
        ts_start[5789] = 1431910950;
        ts_start[5790] = 1431911132;
        ts_start[5791] = 1431913194;
        ts_start[5792] = 1431919249;
        ts_start[5793] = 1431920976;
        ts_start[5794] = 1431922794;
        ts_start[5795] = 1431923711;
        ts_start[5796] = 1431929436;
        ts_start[5797] = 1431936558;
        ts_start[7101] = 1454953454;
        ts_start[7150] = 1455148854;
        ts_start[7263] = 1455508792;
        ts_start[7264] = 1455510422;
        ts_start[7265] = 1455510767;
        ts_start[7266] = 1455511155;
        ts_start[7267] = 1455511582;
        ts_start[7268] = 1455512044;
        ts_start[7271] = 1455517256;
        ts_start[7272] = 1455517571;
        ts_start[7273] = 1455517806;
        ts_start[7274] = 1455518173;
        ts_start[7275] = 1455519253;
        ts_start[7277] = 1455520191;
        ts_start[7278] = 1455521568;
        ts_start[7279] = 1455521994;
        ts_start[7281] = 1455524399;
        ts_start[7291] = 1455528723;
        ts_start[7297] = 1455533692;
        ts_start[7360] = 1455949958;
        ts_start[7361] = 1455969661;
        ts_start[7362] = 1455970571;
        ts_start[7365] = 1455972884;
        ts_start[7366] = 1455973731;
        ts_start[7367] = 1455974320;
        ts_start[7368] = 1455982130;
        ts_start[7370] = 1455986909;
        ts_start[7373] = 1455987766;
        ts_start[7374] = 1455989238;
        ts_start[7396] = 1455999961;
        ts_start[7397] = 1456000200;
        ts_start[7399] = 1456001350;
        ts_start[7401] = 1456002733;
        ts_start[7402] = 1456002959;
        ts_start[7403] = 1456003306;
        ts_start[7404] = 1456004001;
        ts_start[7406] = 1456004665;
        ts_start[7407] = 1456005557;
        ts_start[7408] = 1456006936;
        ts_start[7409] = 1456010987;
        ts_start[7411] = 1456012510;
        ts_start[7412] = 1456016289;
        ts_start[7414] = 1456018683;
        ts_start[7415] = 1456026657;
        ts_start[7418] = 1456027722;
        ts_start[7421] = 1456029670;
        ts_start[7427] = 1456039809;
        ts_start[7428] = 1456041429;
        ts_start[7429] = 1456042300;
        ts_start[7443] = 1456053701;
        ts_start[7444] = 1456056792;
        ts_start[7445] = 1456059051;
        ts_start[7446] = 1456061612;
        ts_start[7447] = 1456067232;
        ts_start[7448] = 1456073194;
        ts_start[7451] = 1456075440;
        ts_start[7453] = 1456076126;
        ts_start[7455] = 1456081492;
        ts_start[7456] = 1456081768;
        ts_start[7457] = 1456082278;
        ts_start[7460] = 1456083354;
        ts_start[7461] = 1456083581;
        ts_start[7464] = 1456084772;
        ts_start[7465] = 1456084972;
        ts_start[7467] = 1456085978;
        ts_start[7468] = 1456086353;
        ts_start[7470] = 1456087146;
        ts_start[7474] = 1456088905;
        ts_start[7475] = 1456089125;
        ts_start[7476] = 1456089440;
        ts_start[7477] = 1456093211;
        ts_start[7478] = 1456098982;
        ts_start[7479] = 1456099490;
        ts_start[7481] = 1456103197;
        ts_start[7482] = 1456103483;
        ts_start[7483] = 1456103865;
        ts_start[7485] = 1456105343;
        ts_start[7486] = 1456105767;
        ts_start[7487] = 1456108133;
        ts_start[7488] = 1456109993;
        ts_start[7489] = 1456112361;
        ts_start[7490] = 1456112584;
        ts_start[7491] = 1456114760;
        ts_start[7492] = 1456114998;
        ts_start[7493] = 1456116899;
        ts_start[7494] = 1456117105;
        ts_start[7495] = 1456120762;
        ts_start[7496] = 1456128846;
        ts_start[7567] = 1456515745;
        ts_start[7579] = 1456526207;
        ts_start[7580] = 1456527004;
        ts_start[7581] = 1456528493;
        ts_start[7589] = 1456592905;
        ts_start[7595] = 1456594992;
        ts_start[7606] = 1456608506;
        ts_start[7607] = 1456609429;
        ts_start[7609] = 1456610461;
        ts_start[7611] = 1456614171;
        ts_start[7612] = 1456615136;
        ts_start[7613] = 1456615554;
        ts_start[7614] = 1456615951;
        ts_start[7616] = 1456624112;
        ts_start[7622] = 1456646964;
        ts_start[7624] = 1456656137;
        ts_start[7628] = 1456672829;
        ts_start[7629] = 1456678135;
        ts_start[7630] = 1456683358;
        ts_start[7631] = 1456688710;
        ts_start[7634] = 1456708998;
        ts_start[7636] = 1456710885;
        ts_start[7637] = 1456722195;
        ts_start[7638] = 1456724456;
        ts_start[7639] = 1456724742;
        ts_start[7641] = 1456726217;
        ts_start[7642] = 1456727243;
        ts_start[7643] = 1456727928;
        ts_start[7644] = 1456728220;
        ts_start[7646] = 1456738877;
        ts_start[7649] = 1456740945;
        ts_start[7652] = 1456742238;
        ts_start[7653] = 1456742432;
        ts_start[7667] = 1456789900;
        ts_start[7669] = 1456790358;
        ts_start[7670] = 1456790499;
        ts_start[7671] = 1456795236;
        ts_start[7779] = 1457153871;
        ts_start[7780] = 1457162724;
        ts_start[7781] = 1457171011;
        ts_start[7782] = 1457181097;
        ts_start[7783] = 1457187461;
        ts_start[7784] = 1457188271;
        ts_start[7785] = 1457190672;
        ts_start[7786] = 1457192984;
        ts_start[7794] = 1457255185;
        ts_start[7795] = 1457257969;
        ts_start[7796] = 1457264901;
        ts_start[7797] = 1457273318;
        ts_start[7798] = 1457274993;
        ts_start[7799] = 1457282999;
        ts_start[7800] = 1457290149;
        ts_start[7801] = 1457297762;
        ts_start[7802] = 1457308873;
        ts_start[7803] = 1457314710;
        ts_start[7804] = 1457322931;
        ts_start[7805] = 1457331200;
        ts_start[7806] = 1457339218;
        ts_start[7807] = 1457339467;
        ts_start[7808] = 1457345498;
        ts_start[7809] = 1457347027;
        ts_start[7946] = 1460209497;
        ts_start[7947] = 1460210764;
        ts_start[7948] = 1460217449;
        ts_start[7949] = 1460225339;
        ts_start[7953] = 1460237473;
        ts_start[7961] = 1460242568;
        ts_start[7962] = 1460245027;
        ts_start[7963] = 1460247048;
        ts_start[7964] = 1460253941;
        ts_start[7965] = 1460261866;
        ts_start[7966] = 1460265770;
        ts_start[7967] = 1460273433;
        ts_start[7968] = 1460280593;
        ts_start[7969] = 1460289043;
        ts_start[7970] = 1460291627;
        ts_start[7971] = 1460300774;
        ts_start[7972] = 1460301629;
        ts_start[7973] = 1460309105;
        ts_start[7975] = 1460310370;
        ts_start[7976] = 1460311551;
        ts_start[7982] = 1460318312;
        ts_start[7983] = 1460320923;
        ts_start[7984] = 1460328303;
        ts_start[7985] = 1460335512;
        ts_start[7986] = 1460343845;
        ts_start[7987] = 1460351451;
        ts_start[7988] = 1460358526;
        ts_start[7989] = 1460366324;
        ts_start[8018] = 1460743793;
        ts_start[8023] = 1460749468;
        ts_start[8024] = 1460792937;
        ts_start[8025] = 1460796016;
        ts_start[8026] = 1460804710;
        ts_start[8027] = 1460811724;
        ts_start[8028] = 1460819012;
        ts_start[8029] = 1460827507;
        ts_start[8030] = 1460835157;
        ts_start[8031] = 1460840901;
        ts_start[8032] = 1460845018;
        ts_start[8033] = 1460845274;
        ts_start[8034] = 1460845969;
        ts_start[8035] = 1460846799;
        ts_start[8036] = 1460847318;
        ts_start[8038] = 1460847627;
        ts_start[8039] = 1460853657;
        ts_start[8040] = 1460860919;
        ts_start[8041] = 1460867852;
        ts_start[8043] = 1460871761;
        ts_start[8044] = 1460878771;
        ts_start[8045] = 1460885664;
        ts_start[8046] = 1460896347;
        ts_start[8047] = 1460904094;
        ts_start[8048] = 1460912108;
        ts_start[8049] = 1460919583;
        ts_start[8050] = 1460931725;
        ts_start[8051] = 1460932346;
        ts_start[8052] = 1460935893;
        ts_start[8053] = 1460941349;
        ts_start[8054] = 1460941822;
        ts_start[8055] = 1460946133;
        ts_start[8056] = 1460951324;
        ts_start[8057] = 1460952094;
        ts_start[8058] = 1460959152;
        ts_start[8059] = 1460966250;
        ts_start[8063] = 1461336746;
        ts_start[8064] = 1461345134;
        ts_start[8066] = 1461345696;
        ts_start[8068] = 1461349403;
        ts_start[8069] = 1461377650;
        ts_start[8071] = 1461379622;
        ts_start[8072] = 1461382255;
        ts_start[8073] = 1461390672;
        ts_start[8074] = 1461399439;
        ts_start[8075] = 1461407690;
        ts_start[8076] = 1461411411;
        ts_start[8077] = 1461411824;
        ts_start[8078] = 1461416796;
        ts_start[8079] = 1461418388;
        ts_start[8081] = 1461425555;
        ts_start[8084] = 1461464737;
        ts_start[8085] = 1461465317;
        ts_start[8086] = 1461470449;
        ts_start[8087] = 1461478081;
        ts_start[8088] = 1461487746;
        ts_start[8089] = 1461493938;
        ts_start[8090] = 1461502480;
        ts_start[8091] = 1461505867;
        ts_start[8092] = 1461508304;
        ts_start[8093] = 1461518427;
        ts_start[8094] = 1461522237;
        ts_start[8095] = 1461530511;
        ts_start[8096] = 1461539023;
        ts_start[8097] = 1461547618;
        ts_start[8098] = 1461555245;
        ts_start[8099] = 1461563175;
        ts_start[8100] = 1461573897;
    }
    
    Long tiTimeOffset;
    
    @Override
    public void process(EventHeader event){
        
        double timestamp = event.getTimeStamp()/1e9; //convert to seconds
        if(timestamp < 1.3e9){ //ie, if this file was reconed before the Ti time offset was put in the database
            if(tiTimeOffset == null){
                tiTimeOffset = getTiTimeOffset(event.getRunNumber());
            }
            timestamp += tiTimeOffset/1e9;
        }
        Map<Cluster, Cluster> replacementClusters = new LinkedHashMap();  //I hate doing this, but for some reason SIOCluster is sometimes used, which does not have a "setEnergy" method.  
        for(Cluster c: event.get(Cluster.class, ecalClusterCollection)){
            if(c instanceof BaseCluster){
                BaseCluster cc = (BaseCluster) c;
                //BaseCluster cc = (SIOCluster) c;

                //double energy = cc.getEnergy();
                double energyCorr = correctedEnergy(c, timestamp);


                cc.setEnergy(energyCorr);
            }
            else{
                BaseCluster replacement = new BaseCluster(c);
                replacement.setNeedsPropertyCalculation(false);
                double energyCorr = correctedEnergy(c, timestamp);
                replacement.setEnergy(energyCorr);
                replacementClusters.put(c, replacement);
               
            }

        }

        
        //replace the original clusters with the replacement versions, and any objects that make reference to clusters.  
        int nClusters = event.get(Cluster.class, ecalClusterCollection).size();
        event.get(Cluster.class, ecalClusterCollection).clear();
        for(Cluster c: replacementClusters.values()){
            event.get(Cluster.class, ecalClusterCollection).add(c);
        }
        
        //check if same number of clusters.  
        if(event.get(Cluster.class, ecalClusterCollection).size() != nClusters){
            throw new RuntimeException("I must have done something wrong:  start with " + nClusters + " clusters, end with " + event.get(Cluster.class, ecalClusterCollection).size() + " clusters" );
        }
        
        
        
        Map<ReconstructedParticle, ReconstructedParticle> replacementParticles = new LinkedHashMap();


        for(ReconstructedParticle particle : event.get(ReconstructedParticle.class, "FinalStateParticles")){


            Cluster replacementCluster = particle.getClusters().size() != 0 ? replacementClusters.get(particle.getClusters().get(0)) : null;

            
            BaseReconstructedParticle replacementParticle = new BaseReconstructedParticle();
            replacementParticle.set4Vector(new BasicHepLorentzVector(replacementCluster != null ? replacementCluster.getEnergy() : 0,particle.getMomentum().v()));
            if(replacementCluster != null)
                replacementParticle.addCluster(replacementCluster);
            replacementParticle.getTracks().addAll(particle.getTracks());
            replacementParticle.setCharge(particle.getCharge());
            replacementParticle.setParticleIdUsed(particle.getParticleIDUsed());
            replacementParticle.setGoodnessOfPid(particle.getGoodnessOfPID());
            replacementParticle.setMass(particle.getMass());
            replacementParticle.setReferencePoint(particle.getReferencePoint());
            replacementParticle.setStartVertex(particle.getStartVertex());
            replacementParticle.setType(particle.getType());
            replacementParticles.put(particle, replacementParticle);

        }
        
        int nFSParticles = event.get(ReconstructedParticle.class, "FinalStateParticles").size();
        event.get(ReconstructedParticle.class, "FinalStateParticles").clear();
        event.get(ReconstructedParticle.class, "FinalStateParticles").addAll(replacementParticles.values());
        if(event.get(ReconstructedParticle.class, "FinalStateParticles").size() != nFSParticles){
            throw new RuntimeException("I must have done something wrong here:  start with " + nFSParticles + " particles, end with " 
                    + event.get(ReconstructedParticle.class, ecalClusterCollection).size() + " particles" );
        }
            
        //now for the v0 and moller collections
        
        outerLoop : for(List<ReconstructedParticle> particleList : event.get(ReconstructedParticle.class)){
            Map<ReconstructedParticle, ReconstructedParticle> replacementVertexParticles = new LinkedHashMap();  // a map containing only the replacements for this type of vertex
            for(ReconstructedParticle particle : particleList){
                if(particle.getParticles().size() != 2) 
                    continue outerLoop;
                ReconstructedParticle p1 = replacementParticles.get(particle.getParticles().get(0));
                ReconstructedParticle p2 = replacementParticles.get(particle.getParticles().get(1));
                if(p1 == null)  //happens if one of the particles didn't need to be replaced.  
                    p1 = particle.getParticles().get(0);
                if(p2 == null) 
                    p2 = particle.getParticles().get(1);

                BaseReconstructedParticle replacementParticle = new BaseReconstructedParticle();
                double esum = p1.getEnergy() + p2.getEnergy();
                Hep3Vector v = particle.getMomentum();
                replacementParticle.set4Vector(new BasicHepLorentzVector(esum,v));
                replacementParticle.addParticle(p1);
                replacementParticle.addParticle(p2);
                replacementParticle.setCharge(particle.getCharge());
 //               replacementParticle.setParticleIdUsed(particle.getParticleIDUsed());
                replacementParticle.setGoodnessOfPid(particle.getGoodnessOfPID());
                replacementParticle.setMass(particle.getMass());
                replacementParticle.setReferencePoint(particle.getReferencePoint());
                replacementParticle.setStartVertex(particle.getStartVertex());
                replacementParticle.setType(particle.getType());
                replacementVertexParticles.put(particle, replacementParticle);

            }
            particleList.removeAll(replacementVertexParticles.keySet());
            particleList.addAll(replacementVertexParticles.values());
            replacementParticles.putAll(replacementVertexParticles);
        }

        //now for the vertices formed from the moller and v0 particle objects. 
        
        Map<Vertex, Vertex> allReplacementVertices = new LinkedHashMap();
        for(List<Vertex> vertexList : event.get(Vertex.class)){
            List<Vertex> replacementVertices = new ArrayList();
            for(Vertex v : vertexList){
                BaseVertex bv = new BaseVertex(v.isPrimary(), v.getAlgorithmType(), v.getChi2(), v.getProbability(), v.getCovMatrix(), v.getPosition(), 
                        replacementParticles.get(v.getAssociatedParticle()));
                replacementVertices.add(bv);
                allReplacementVertices.put(v, bv);
            }
            vertexList.clear();
            vertexList.addAll(replacementVertices);
        }
        for(List<ReconstructedParticle> particleList : event.get(ReconstructedParticle.class)){
            for(ReconstructedParticle particle : particleList){
                ((BaseReconstructedParticle)particle).setStartVertex(allReplacementVertices.get(particle.getStartVertex()));
            }
        }


    }
    Logger LOGGER = Logger.getLogger(TimeDependentGainTweakDriver.class.getPackage().getName());


    private final long timestampCycle = 24 * 6 * 35;
    private Long getTiTimeOffset(int run) {
        Long currentTiTimeOffset = null;
        RunManager runManager = RunManager.getRunManager();
       // runManager.setRun(runNumber);
        if (runManager.getRun() != null) {
            if (runManager.runExists()) {
                currentTiTimeOffset = runManager.getRunSummary().getTiTimeOffset();
                tiTimeOffset = (currentTiTimeOffset / timestampCycle) * timestampCycle;
                LOGGER.info("TI time offset set to " + currentTiTimeOffset + " for run "
                        + run + " from database");
            } else {
                LOGGER.warning("Run " + run 
                        + " does not exist in the run database.");
            }
        } else {
            LOGGER.info("Run manager is not initialized; TI time offset not available.");
        }
        /* Make sure connection is closed immediately. --JM */
        try {
            LOGGER.info("Closing run manager db connection ...");
            RunManager.getRunManager().closeConnection();
            LOGGER.info("Run manager db connection was closed.");
        } catch (Exception e) {
            e.printStackTrace();
        }
       return tiTimeOffset;
    }



    private double correctedEnergy(Cluster c, double timestamp) {
        double totHitEnergy = 0;
        double corrTotHitEnergy = 0;
        for(CalorimeterHit hit : c.getCalorimeterHits()){
            double energy = hit.getCorrectedEnergy();
            totHitEnergy += energy;
            corrTotHitEnergy += energy * getCorrFactor(hit.getIdentifierFieldValue("ix"), hit.getIdentifierFieldValue("iy"), timestamp);
        }
        return c.getEnergy()*corrTotHitEnergy/totHitEnergy;
    }



    abstract double getCorrFactor(int ix, int iy, double timestamp);


}
