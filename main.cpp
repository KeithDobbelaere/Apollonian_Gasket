#include <SFML/Graphics.hpp>
#include <cmath>
#include <complex>
#include <vector>
#include <array>
#include <iostream>

#define M_PI 3.14159265358979323846
static const double epsilon = 0.1; 

class Circle : public sf::Drawable
{
public:
    Circle(double curvature, const std::complex<double> &center, const sf::Color &color = sf::Color::White)
    : m_curvature(curvature), m_radius(abs(1.0/curvature)), m_center(center), m_circleShape((float)m_radius)
    {
        const float smoothnessFactor = 4.f;
        auto pointCount = (size_t)std::max(20.f, sqrt((float)m_radius) * smoothnessFactor);
        m_circleShape.setPointCount(pointCount);
        m_circleShape.setOrigin({m_circleShape.getRadius(), m_circleShape.getRadius()});
        m_circleShape.setPosition({(float)center.real(), (float)center.imag()});
        m_circleShape.setFillColor(sf::Color::Transparent);
        m_circleShape.setOutlineColor(color);
        m_circleShape.setOutlineThickness(-1);
    }

    double getRadius() const {
        return m_radius;
    }

    double getCurvature() const {
        return m_curvature;
    }

    const std::complex<double> &getCenter() const {
        return m_center;
    }

    double distance(const Circle &other) const {
        double real = m_center.real() - other.getCenter().real();
        double imag = m_center.imag() - other.getCenter().imag();
        return abs(sqrt(real * real + imag * imag));
    }

private:
    void draw(sf::RenderTarget &target, const sf::RenderStates &states) const override {
        target.draw(m_circleShape, states);
    }

    double m_curvature;
    double m_radius;
    std::complex<double> m_center;
    sf::CircleShape m_circleShape;
};

class ApollonianGasket : public sf::Drawable, public sf::Transformable
{
public:
    ApollonianGasket() {
        // Randomize seed
        const auto seed = time(NULL);
        #ifndef NDEBUG
        std::cout << "Seed: " << seed << std::endl;
        #endif
        srand(seed);

        initialize();
    }

    void initialize() {
        const double C1_RADIUS = 400.0;
        // Create first circle
        auto c1 = Circle(-1.0/C1_RADIUS, {C1_RADIUS, C1_RADIUS});
        // Get random radius between 20 and c1's radius / 2
        double r2 = (double)(rand() % 200 + 20);
        // Get random position for c2 tangent to c1
        auto v = random2DVector();
        vectorSetMagnitude(v, c1.getRadius() - r2);
        auto c2 = Circle(1.0/r2, {C1_RADIUS + v.x, C1_RADIUS + v.y});
        // Get size and position for c3 tangent to c1 and c2
        double r3 = v.length();
        v = v.rotatedBy(sf::degrees(180.0));
        vectorSetMagnitude(v, c1.getRadius() - r3);
        auto c3 = Circle(1.0/r3, {C1_RADIUS + v.x, C1_RADIUS + v.y}, sf::Color::Yellow);

        m_circles = {c1, c2, c3};
        m_queue = {{c1, c2, c3}};
    }

    void findNextCircles() {
        std::vector<std::tuple<Circle, Circle, Circle>> nextQueue;
        int limit = 10000;
        for (auto &triplet : m_queue)
        {
            auto &[c1, c2, c3] = triplet;

            const auto &k4 = descartes(c1, c2, c3);
            const auto &newCircles = complexDescartes(c1, c2, c3, k4);
            for (const auto &circle : newCircles)
            {
                if (validate(circle, c1, c2, c3))
                {
                    m_circles.push_back(circle);
                    nextQueue.emplace_back(c1, c2, m_circles.back());
                    nextQueue.emplace_back(c1, c3, m_circles.back());
                    nextQueue.emplace_back(c2, c3, m_circles.back());
                }
            }
            if (limit-- < 0)
                break;
        }
        m_queue = nextQueue;
    }

private:
    std::array<Circle, 4> complexDescartes(const Circle &c1, const Circle &c2, const Circle &c3, const std::pair<double, double> &k4)
    {
        const double &k1 = c1.getCurvature();
        const double &k2 = c2.getCurvature();
        const double &k3 = c3.getCurvature();

        const auto &z1 = c1.getCenter();
        const auto &z2 = c2.getCenter();
        const auto &z3 = c3.getCenter();

        auto zk1 = z1 * k1;
        auto zk2 = z2 * k2;
        auto zk3 = z3 * k3;
        auto sum = zk1 + zk2 + zk3;

        auto root = zk1 * zk2 + zk2 * zk3 + zk1 * zk3;
        root = sqrt(root) * 2.0;

        auto z4_1 = (sum + root) * (1.0 / k4.first);
        auto z4_2 = (sum - root) * (1.0 / k4.first);
        auto z4_3 = (sum + root) * (1.0 / k4.second);
        auto z4_4 = (sum - root) * (1.0 / k4.second);
        return {
            Circle(k4.first,  z4_1),
            Circle(k4.first,  z4_2),
            Circle(k4.second, z4_3),
            Circle(k4.second, z4_4)}; 
    }

    std::pair<double, double> descartes(const Circle &c1, const Circle &c2, const Circle &c3)
    {
        const double &k1 = c1.getCurvature();
        const double &k2 = c2.getCurvature();
        const double &k3 = c3.getCurvature();

        double sum = k1 + k2 + k3;
        double product = abs(k1 * k2 + k2 * k3 + k1 * k3);
        double root = 2.0 * sqrt(product);

        return {sum + root, sum - root};
    }

    bool isTangent(const Circle &c1, const Circle &c2)
    {
        const double &distance = c1.distance(c2);
        const double &r1 = c1.getRadius();
        const double &r2 = c2.getRadius();

        bool a = abs(distance - (r1 + r2)) < epsilon;
        bool b = abs(distance - abs(r1 - r2)) < epsilon;
        #ifndef NDEBUG 
        if (!(a || b)) {
            // Print rejected circle info
            Circle badCircle(c1.getCurvature(), c1.getCenter(), sf::Color::Red);
            m_circles.emplace_back(badCircle);
            std::cerr << "Rejected circle - position: " << c1.getCenter() << ", radius: " << c1.getRadius() << "\n";
        }
        #endif
        return (a || b);
    }

    bool validate(const Circle &c4, const Circle &c1, const Circle &c2, const Circle &c3)
    {
        // Eliminate circles that are too small
        if (c4.getRadius() < 1.0)
            return false;
        // Eliminate duplicates
        for (const auto &other : m_circles) {
            const double &distance = c4.distance(other);
            double rDiff = abs(c4.getRadius() - other.getRadius());
            if (distance < epsilon && rDiff < epsilon)
                return false;
        }
        // Eliminate circles that are not tangent to the other three
        if (!isTangent(c4, c1) || !isTangent(c4, c2) || !isTangent(c4, c3))
            return false;

        return true;
    }

    sf::Vector2<double> random2DVector()
    {
        double angle = (double)(rand()) / ((double)(RAND_MAX / (2.0 * M_PI)));
        const double length = 1.0;
        return {cos(angle) * length, sin(angle) * length};
    }

    void vectorSetMagnitude(sf::Vector2<double>& v, double length)
    {
        double angle = atan2(v.y, v.x);
        v = sf::Vector2<double>(cos(angle) * length, sin(angle) * length);
    }

    void draw(sf::RenderTarget &target, const sf::RenderStates& states) const override
    {
        const auto transform = getTransform() * states.transform;

        for (const auto &circle : m_circles)
        {
            target.draw(circle, transform);
        }
    }

    std::vector<Circle> m_circles;
    std::vector<std::tuple<Circle, Circle, Circle>> m_queue;
};

void mousePressed(const sf::Event &event, ApollonianGasket &apg)
{
    if (event.mouseButton.button == sf::Mouse::Button::Left)
        apg.findNextCircles();
    if (event.mouseButton.button == sf::Mouse::Button::Right)
        apg.initialize();
}

constexpr unsigned int DEFAULT_WIDTH = 800;
constexpr unsigned int DEFAULT_HEIGHT = 800;
constexpr float ASPECT_RATIO = float(DEFAULT_WIDTH) / float(DEFAULT_HEIGHT);
void resizeWindow(const sf::Event &event, sf::RenderWindow &window)
 {
    // Maintain aspect ratio
    unsigned newWidth = event.size.width;
    unsigned newHeight = event.size.height;
    if (newWidth / ASPECT_RATIO < newHeight)
        newHeight = newWidth / ASPECT_RATIO;
    else
        newWidth = newHeight * ASPECT_RATIO; 
    window.setSize({newWidth, newHeight});
}

int main()
{
    sf::RenderWindow window(sf::VideoMode({DEFAULT_WIDTH, DEFAULT_HEIGHT}), "Apolloninan Gasket");
    window.setFramerateLimit(60);

    ApollonianGasket apg;
    apg.setPosition({DEFAULT_WIDTH / 2.f - 400.f, DEFAULT_HEIGHT / 2.f - 400.f});

    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
            if (event.type == sf::Event::Resized)
                resizeWindow(event, window);
            if (window.hasFocus() && event.type == sf::Event::MouseButtonPressed)
                mousePressed(event, apg);
        }
        //ap.findNextCircles();
        window.clear();
        window.draw(apg);
        window.display();
    }

    return 0;
}